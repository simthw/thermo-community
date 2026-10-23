module Analysis

using CSV
using DataFrames
using Dates
using DifferentialEquations
using Distributed
using Distributions

using JLD2
using Printf
using Random
using SHA
using Statistics
using TOML

include("constants.jl")
include("assembly_cook/AssemblyCook.jl")
using .AssemblyCook

export anmbrRunConfig, run_anmbr_analysis, run_anmbr_analysis_parallel, ensure_output_tree, DEFAULT_OUTPUT_ROOT
export MODEL_REVISION, OUTPUT_SCHEMA_VERSION
export model_input_hashes
export initialise_anmbr, initialise_adm1_like, full_simulate
export run_single_simulation  # Export single simulation function for distributed use

# Export from assembly_cook
export Reaction, make_Reaction
export MicrobeP, make_MicrobeP
export FullParameters, make_FullParameters, extinction
export initialise, initialise_syn_trio, initialise_pol_duo
export full_simulate_syn
export θ, qs
export calculate_reaction_gibbs, analyze_thermodynamic_feasibility
export calculate_actual_ΔG, reaction_entropy_production, total_entropy_production
export get_reaction_energy_sum
export options_titles, initialise_np, initialise_Pb, initialise_ϕH, initialise_sat
export Force, nForce, Jacobian, nJacobian
export quantify_ints, analyse_ints, InteractionResult
export robustness_test, batch_robustness_test, RobustnessResult

const PROJECT_ROOT = rstrip(normpath(joinpath(@__DIR__, "..")), ['\\', '/'])
const DEFAULT_OUTPUT_ROOT = normpath(joinpath(PROJECT_ROOT, "simulations"))
const MODEL_REVISION = "2026-07-22"
const OUTPUT_SCHEMA_VERSION = 3
const MODEL_PROVENANCE_PATHS = vcat(
  [joinpath(@__DIR__, "constants.jl")],
  sort(filter(
    path -> endswith(path, ".jl"),
    readdir(joinpath(@__DIR__, "assembly_cook"); join=true),
  )),
  [
    joinpath(@__DIR__, "data", "reactions.dat"),
    joinpath(@__DIR__, "data", "metabolite_energies.csv"),
  ],
)

function model_input_hashes()
  return Dict(
    replace(relpath(path, PROJECT_ROOT), '\\' => '/') => bytes2hex(SHA.sha256(read(path)))
    for path in MODEL_PROVENANCE_PATHS
  )
end

function project_output_path(path::AbstractString)
  isabspath(path) && return normpath(path)
  normalized = normpath(String(path))
  project_name = basename(PROJECT_ROOT)
  normalized_unix = replace(normalized, '\\' => '/')
  if normalized_unix == project_name || startswith(normalized_unix, project_name * "/")
    return normpath(joinpath(dirname(PROJECT_ROOT), normalized))
  end
  return normpath(joinpath(PROJECT_ROOT, normalized))
end

"""
    anmbrRunConfig(; kwargs...)

Container for the configuration required to reproduce the Assembly simulations
with the `initialise_anmbr` parameter generator.

Keyword arguments:

- `repeats::Int = 5`: number of independent parameter sets to generate and simulate.
- `kc::Float64 = 10.0`: average catalytic rate constant passed to `initialise_anmbr`.
- `KS::Float64 = (1 / 4) * 5.5e-3`: average half-saturation constant passed to `initialise_anmbr`.
- `kr::Float64 = 10.0`: average reversibility constant passed to `initialise_anmbr`.
- `Tmax::Float64 = 2.0e6`: maximum integration time (seconds) for `full_simulate`.
- `initial_population::Float64 = 1e5`: starting cell count for each strain.
- `initial_energy::Float64 = 1.0e5`: initial ATP-equivalent reserves for each strain.
- `initial_fraction::Float64 = 0.128`: initial ribosome fraction for each strain.
- `glucose_mM::Float64 = 5.0`: initial glucose concentration in mM.
- `alkalinity_mgCaCO3_L::Float64 = 1200.0`: alkalinity as mg CaCO3/L (contributes HCO3-).
- `sulfate_mg_L::Float64 = 30.0`: sulfate concentration in mg/L (baseline scenario).
- `seed::Int = 1`: base seed controlling the random stream used by `initialise_anmbr`.
- `output_root::AbstractString = DEFAULT_OUTPUT_ROOT`: directory where run artefacts are stored.
- `save_trajectories::Bool = true`: toggle for persisting full time series to disk.
- `thermodynamic_feedback::Bool = true`: include theta in the reaction rate law.
"""
struct anmbrRunConfig
  repeats::Int
  kc::Float64
  KS::Float64
  kr::Float64
  Tmax::Float64
  initial_population::Float64
  initial_energy::Float64
  initial_fraction::Float64
  glucose_mM::Float64
  alkalinity_mgCaCO3_L::Float64
  sulfate_mg_L::Float64
  seed::Int
  output_root::String
  save_trajectories::Bool
  thermodynamic_feedback::Bool
end

function anmbrRunConfig(; repeats::Int=5, kc::Float64=10.0,
  KS::Float64=(1 / 4) * 5.5e-3, kr::Float64=10.0, Tmax::Float64=2.0e6,
  initial_population::Float64=1e5, initial_energy::Float64=1.0e5,
  initial_fraction::Float64=0.128, glucose_mM::Float64=5.0,
  alkalinity_mgCaCO3_L::Float64=1200.0, sulfate_mg_L::Float64=30.0,
  seed::Int=1, output_root::AbstractString=DEFAULT_OUTPUT_ROOT,
  save_trajectories::Bool=true, thermodynamic_feedback::Bool=true)
  repeats > 0 || throw(ArgumentError("repeats must be positive"))
  kc > 0 || throw(ArgumentError("kc must be positive"))
  KS >= 0 || throw(ArgumentError("KS must be non-negative"))
  kr > 0 || throw(ArgumentError("kr must be positive"))
  Tmax > 0 || throw(ArgumentError("Tmax must be positive"))
  initial_population >= 0 || throw(ArgumentError("initial_population must be non-negative"))
  initial_energy >= 0 || throw(ArgumentError("initial_energy must be non-negative"))
  0.0 <= initial_fraction <= 1.0 || throw(ArgumentError("initial_fraction must be within [0, 1]"))
  glucose_mM >= 0 || throw(ArgumentError("glucose_mM must be non-negative"))
  alkalinity_mgCaCO3_L >= 0 || throw(ArgumentError("alkalinity must be non-negative"))
  sulfate_mg_L >= 0 || throw(ArgumentError("sulfate must be non-negative"))

  return anmbrRunConfig(
    repeats,
    kc,
    KS,
    kr,
    Tmax,
    initial_population,
    initial_energy,
    initial_fraction,
    glucose_mM,
    alkalinity_mgCaCO3_L,
    sulfate_mg_L,
    seed,
    project_output_path(output_root),
    save_trajectories,
    thermodynamic_feedback,
  )
end

"""
    ensure_output_tree(config)

Create the directory structure used by the anmbr pipeline for a given `config`
and return a named tuple with the realised paths.
"""
function ensure_output_tree(config::anmbrRunConfig)
  run_stamp = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHHMMSSsss") * "Z"
  run_dir = joinpath(config.output_root, run_stamp)
  paras_dir = joinpath(run_dir, "paras")
  traj_dir = joinpath(run_dir, "trajectories")

  mkpath(config.output_root)
  ispath(run_dir) && error("Refusing to reuse existing run directory: $run_dir")
  mkdir(run_dir)
  mkdir(paras_dir)
  if config.save_trajectories
    mkdir(traj_dir)
  end

  return (run_dir=run_dir, paras_dir=paras_dir, traj_dir=traj_dir)
end

function write_config(config::anmbrRunConfig, run_dir::String, seeds::Vector{Int})
  config_table = Dict(
    "repeats" => config.repeats,
    "kc" => config.kc,
    "KS" => config.KS,
    "kr" => config.kr,
    "Tmax" => config.Tmax,
    "initial_population" => config.initial_population,
    "initial_energy" => config.initial_energy,
    "initial_fraction" => config.initial_fraction,
    "glucose_mM" => config.glucose_mM,
    "alkalinity_mgCaCO3_L" => config.alkalinity_mgCaCO3_L,
    "sulfate_mg_L" => config.sulfate_mg_L,
    "seed" => config.seed,
    "seed_sequence" => seeds,
    "output_root" => config.output_root,
    "save_trajectories" => config.save_trajectories,
    "thermodynamic_feedback" => config.thermodynamic_feedback,
    "schema_version" => OUTPUT_SCHEMA_VERSION,
    "model_revision" => MODEL_REVISION,
    "input_sha256" => model_input_hashes(),
    "julia_version" => string(VERSION),
    "generated_at_utc" => string(now(UTC)),
  )
  open(joinpath(run_dir, "config.toml"), "w") do io
    TOML.print(io, config_table)
  end
  return nothing
end

function write_completion_manifest(run_dir::String, summary::DataFrame)
  successful = :solver_success in propertynames(summary) ? all(Bool.(summary.solver_success)) : false
  reached_tmax = :reached_tmax in propertynames(summary) ? all(Bool.(summary.reached_tmax)) : false
  manifest = Dict(
    "schema_version" => OUTPUT_SCHEMA_VERSION,
    "model_revision" => MODEL_REVISION,
    "input_sha256" => model_input_hashes(),
    "completed_at_utc" => string(now(UTC)),
    "runs" => nrow(summary),
    "all_solver_success" => successful,
    "all_reached_tmax" => reached_tmax,
  )
  open(joinpath(run_dir, "COMPLETED.toml"), "w") do io
    TOML.print(io, manifest)
  end
  return nothing
end

"""
    run_anmbr_analysis([config::anmbrRunConfig])

Execute a batch of simulations using `initialise_anmbr` to reproduce the
ANMBR pipeline with the constrained parameter set.

Returns the path to the directory containing the generated artefacts.
"""
function run_anmbr_analysis(config::anmbrRunConfig=anmbrRunConfig())
  paths = ensure_output_tree(config)
  seeds = config.seed .+ (0:(config.repeats-1))
  write_config(config, paths.run_dir, collect(seeds))

  summary_rows = NamedTuple[]
  for (idx, seed) in enumerate(seeds)
    println("⟹ Starting anmbr simulation $(idx)/$(config.repeats) with seed = $(seed)")
    Random.seed!(seed)

    κ_sulfate = AssemblyCook.DEFAULT_d * config.sulfate_mg_L / 96000.0
    κ_HCO3 = AssemblyCook.DEFAULT_d * config.alkalinity_mgCaCO3_L / 100000.0
    ps = initialise_anmbr(config.kc, config.KS, config.kr;
      κ_sulfate=κ_sulfate, κ_HCO3=κ_HCO3, verbose=false)
    N = ps.N
    M = ps.M

    pop0 = fill(config.initial_population, N)

    # Initialize metabolite concentrations based on actual anaerobic digestion conditions
    # Metabolite indices:
    #   1: Glucose, 2: HCO3-, 3: Acetate, 4: CO2, 5: H2O
    #   6: SO42-, 7: HS-, 8: H+, 9: Lactate, 10: Butyrate
    #   11: Propionate, 12: Valerate, 13: H2, 14: CH4, 15: CO32-
    # Initialized to zero unless set below (Acetate/HS-/Lactate/Butyrate/
    # Propionate/Valerate/H2/CH4 are all produced by fermentation/reduction).
    conc0 = zeros(M)
    conc0[1] = config.glucose_mM / 1000.0                 # Glucose (M)
    conc0[2] = config.alkalinity_mgCaCO3_L / 100000.0     # HCO3- from alkalinity (M)
    conc0[4] = 1.5e-5                                     # CO2: atmospheric equilibrium
    conc0[5] = 1.0                                        # Water: pure-liquid activity ≡ 1 (buffered; enters Q as 1)
    conc0[6] = config.sulfate_mg_L / 96000.0              # SO42- (M) from mg/L
    conc0[8] = 1e-7                                       # H+: pH 7
    conc0[15] = config.alkalinity_mgCaCO3_L / 1000000.0   # CO32-: ~1/10 of HCO3-

    energy0 = fill(config.initial_energy, N)
    phi0 = fill(config.initial_fraction, N)

    t_start = time()
    C, T, solver_metadata = full_simulate(
      ps,
      config.Tmax,
      pop0,
      conc0,
      energy0,
      phi0;
      thermodynamic_feedback=config.thermodynamic_feedback,
      return_metadata=true,
    )
    elapsed = time() - t_start
    solver_metadata.successful || @warn "Solver did not complete successfully" seed solver_metadata

    final_pop = C[end, 1:N]
    final_conc = C[end, (N+1):(N+M)]
    final_energy = C[end, (N+M+1):(N+M+N)]
    final_phi = C[end, (N+M+N+1):(N+M+2N)]

    survivor_mask = final_pop .> 1e-9
    survivors = count(identity, survivor_mask)
    extinct = N - survivors
    biomass = sum(final_pop)
    mean_energy = mean(final_energy)
    mean_fraction = mean(final_phi)
    t_final = T[end]

    param_path = joinpath(paths.paras_dir, @sprintf("anmbr_parameters_run%03d.jld2", idx))
    schema_version = OUTPUT_SCHEMA_VERSION
    model_revision = MODEL_REVISION
    input_sha256 = model_input_hashes()
    thermodynamic_feedback = config.thermodynamic_feedback
    JLD2.@save param_path ps seed pop0 conc0 energy0 phi0 schema_version model_revision input_sha256 thermodynamic_feedback solver_metadata

    if config.save_trajectories
      traj_path = joinpath(paths.traj_dir, @sprintf("anmbr_trajectory_run%03d.jld2", idx))
      JLD2.@save traj_path C T seed final_pop final_conc final_energy final_phi schema_version model_revision input_sha256 thermodynamic_feedback solver_metadata
    end

    push!(summary_rows, (
      run=idx,
      seed=seed,
      survivors=survivors,
      extinct=extinct,
      final_time=t_final,
      runtime_seconds=elapsed,
      biomass=biomass,
      mean_energy=mean_energy,
      mean_fraction=mean_fraction,
      CH4_M=final_conc[14],
      SO4_M=final_conc[6],
      HS_M=final_conc[7],
      solver_retcode=solver_metadata.retcode,
      solver_success=solver_metadata.successful,
      reached_tmax=solver_metadata.reached_tmax,
    ))
  end

  summary_df = DataFrame(summary_rows)
  CSV.write(joinpath(paths.run_dir, "summary.csv"), summary_df)
  write_completion_manifest(paths.run_dir, summary_df)

  println("✓ anmbr analysis finished. Outputs archived in: $(paths.run_dir)")
  return paths.run_dir
end

"""
    run_single_simulation(idx, seed, config, paras_dir, traj_dir)

Execute a single simulation (for serial or parallel calls).

Returns a named tuple containing simulation results.
"""
function run_single_simulation(
  idx::Int,
  seed::Int,
  config::anmbrRunConfig,
  paras_dir::String,
  traj_dir::String
)
  # Set independent random seed
  Random.seed!(seed)

  κ_sulfate = AssemblyCook.DEFAULT_d * config.sulfate_mg_L / 96000.0
  κ_HCO3 = AssemblyCook.DEFAULT_d * config.alkalinity_mgCaCO3_L / 100000.0
  ps = initialise_anmbr(config.kc, config.KS, config.kr;
    κ_sulfate=κ_sulfate, κ_HCO3=κ_HCO3, verbose=false)
  N = ps.N
  M = ps.M

  pop0 = fill(config.initial_population, N)

  # Initialize metabolite concentrations
  conc0 = zeros(M)
  conc0[1] = config.glucose_mM / 1000.0
  conc0[2] = config.alkalinity_mgCaCO3_L / 100000.0
  conc0[4] = 1.5e-5
  conc0[5] = 1.0                                          # Water: pure-liquid activity ≡ 1 (buffered; enters Q as 1)
  conc0[6] = config.sulfate_mg_L / 96000.0
  conc0[8] = 1e-7
  conc0[15] = config.alkalinity_mgCaCO3_L / 1000000.0

  energy0 = fill(config.initial_energy, N)
  phi0 = fill(config.initial_fraction, N)

  t_start = time()
  C, T, solver_metadata = full_simulate(
    ps,
    config.Tmax,
    pop0,
    conc0,
    energy0,
    phi0;
    verbose=false,
    thermodynamic_feedback=config.thermodynamic_feedback,
    return_metadata=true,
  )
  elapsed = time() - t_start
  solver_metadata.successful || @warn "Solver did not complete successfully" seed solver_metadata

  final_pop = C[end, 1:N]
  final_conc = C[end, (N+1):(N+M)]
  final_energy = C[end, (N+M+1):(N+M+N)]
  final_phi = C[end, (N+M+N+1):(N+M+2N)]

  survivor_mask = final_pop .> 1e-9
  survivors = count(identity, survivor_mask)
  extinct = N - survivors
  biomass = sum(final_pop)
  mean_energy = mean(final_energy)
  mean_fraction = mean(final_phi)
  t_final = T[end]

  # Save parameters
  param_path = joinpath(paras_dir, @sprintf("anmbr_parameters_run%03d.jld2", idx))
  schema_version = OUTPUT_SCHEMA_VERSION
  model_revision = MODEL_REVISION
  input_sha256 = model_input_hashes()
  thermodynamic_feedback = config.thermodynamic_feedback
  JLD2.@save param_path ps seed pop0 conc0 energy0 phi0 schema_version model_revision input_sha256 thermodynamic_feedback solver_metadata

  # Save trajectories (if enabled)
  if config.save_trajectories
    traj_path = joinpath(traj_dir, @sprintf("anmbr_trajectory_run%03d.jld2", idx))
    JLD2.@save traj_path C T seed final_pop final_conc final_energy final_phi schema_version model_revision input_sha256 thermodynamic_feedback solver_metadata
  end

  return (
    run=idx,
    seed=seed,
    survivors=survivors,
    extinct=extinct,
    final_time=t_final,
    runtime_seconds=elapsed,
    biomass=biomass,
    mean_energy=mean_energy,
    mean_fraction=mean_fraction,
    CH4_M=final_conc[14],
    SO4_M=final_conc[6],
    HS_M=final_conc[7],
    solver_retcode=solver_metadata.retcode,
    solver_success=solver_metadata.successful,
    reached_tmax=solver_metadata.reached_tmax,
  )
end

"""
    run_anmbr_analysis_parallel([config::anmbrRunConfig]; nworkers::Int=nworkers())

Run simulations in parallel using Distributed.pmap for multi-process parallelism.

# Parameters
- `config`: Configuration struct
- `nworkers`: Number of worker processes (default: all available workers)

# Example
```julia
using Distributed
addprocs(4)  # Add 4 worker processes
@everywhere using Analysis
run_anmbr_analysis_parallel()
```

Return the result directory path.
"""
function run_anmbr_analysis_parallel(config::anmbrRunConfig=anmbrRunConfig(); nworkers::Int=Distributed.nworkers())
  # Check number of worker processes
  actual_workers = Distributed.nworkers()
  if actual_workers < 2
    @warn "Only $(actual_workers) workers available, recommend using addprocs() to add more processes for parallel speedup."
  end

  paths = ensure_output_tree(config)
  seeds = config.seed .+ (0:(config.repeats-1))
  write_config(config, paths.run_dir, collect(seeds))

  println("⟹ Starting parallel simulation: $(config.repeats) runs, $(actual_workers) workers")

  # Execute in parallel using pmap
  # Note: Worker processes need to have Analysis module loaded
  results = Distributed.pmap(
    (args) -> run_single_simulation(args[1], args[2], config, paths.paras_dir, paths.traj_dir),
    [(idx, seed) for (idx, seed) in enumerate(seeds)];
    batch_size=max(1, div(config.repeats, actual_workers))
  )

  # Collect results and write summary
  summary_df = DataFrame(results)
  CSV.write(joinpath(paths.run_dir, "summary.csv"), summary_df)
  write_completion_manifest(paths.run_dir, summary_df)

  println("✓ Parallel simulation complete. Output directory: $(paths.run_dir)")
  return paths.run_dir
end

end # module

#!/usr/bin/env julia
# Deterministic five-guild baseline — 5-strain SO4 sweep.
# Each guild collapsed to ONE deterministic strain carrying ALL guild reactions
# (central kc/KS, midpoint eta, uniform phiP). It retains the metabolite system
# and theta feedback but changes several factors together; it is not ADM1 and
# not a diversity-only intervention.
# Per-strain initial biomass is scaled so each ADM1-like "guild pool" has the
# same total biomass as the corresponding guild's 60-strain ensemble share
# (e.g., methanogen pool starts with 15 × INIT_POP because the ensemble has
# 15 methanogen strains × INIT_POP each).
# Output: Analysis/simulations/ablation_adm1/<timestamp>-SO4_<level>/
#   Same directory layout as `run_anmbr_analysis` so existing post-processing
#   (dissipation_rate.jl) can be re-run on these dirs unchanged.
# Usage:
#   julia --project=Analysis Analysis/src/scripts/run_ablation_adm1.jl

using Analysis
using Dates
using JLD2
using CSV
using DataFrames
using Statistics
using Printf
using TOML

const ABLATION_ROOT = joinpath(Analysis.DEFAULT_OUTPUT_ROOT, "ablation_adm1")
isdir(ABLATION_ROOT) || mkpath(ABLATION_ROOT)

# Numerical defaults mirror anmbrRunConfig
const KC = 10.0
const KS = (1 / 4) * 5.5e-3
const KR = 10.0
const TMAX = 2.0e6
const INIT_POP = 1e5
const INIT_ENERGY = 1.0e5
const INIT_FRAC = 0.128
const GLUCOSE_MM = 5.0
const ALK = 1200.0
const SO4_LEVELS = [10.0, 30.0, 100.0, 300.0, 1000.0]

# Guild sizes in the 60-strain ensemble = round(STRAIN_TYPE_RATIOS × 60).
# Order MUST match `_build_type_to_reactions_map`'s strain_types:
#   1: acidifier, 2: hydrogen_producer, 3: methanogen, 4: sulfate_reducer, 5: acetogen
const GUILD_SIZES = [21, 12, 15, 6, 6]   # sum = 60

ts = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHHMMSSsss") * "Z"
println("ADM1-like baseline sweep — timestamp: $ts")
println("Pop0 per guild = $GUILD_SIZES × $(Int(INIT_POP)) = $(GUILD_SIZES .* INIT_POP)")
println("SO₄²⁻ levels: $SO4_LEVELS mg/L")
println()

created_dirs = String[]
for (i, so4) in enumerate(SO4_LEVELS)
  run_dir = joinpath(ABLATION_ROOT, @sprintf("%s-SO4_%d", ts, Int(so4)))
  paras_dir = joinpath(run_dir, "paras")
  traj_dir = joinpath(run_dir, "trajectories")
  mkpath(paras_dir);
  mkpath(traj_dir)

  κ_sulfate = Analysis.AssemblyCook.DEFAULT_d * so4 / 96000.0
  κ_HCO3 = Analysis.AssemblyCook.DEFAULT_d * ALK / 100000.0
  ps = initialise_adm1_like(KC, KS, KR; κ_sulfate=κ_sulfate, κ_HCO3=κ_HCO3)

  N = ps.N;
  M = ps.M
  @assert N == length(GUILD_SIZES) "ADM1-like N=$N ≠ length(GUILD_SIZES)=$(length(GUILD_SIZES))"

  pop0 = Float64.(GUILD_SIZES) .* INIT_POP
  energy0 = fill(INIT_ENERGY, N)
  phi0 = fill(INIT_FRAC, N)

  conc0 = zeros(M)
  conc0[1] = GLUCOSE_MM / 1000.0
  conc0[2] = ALK / 100000.0
  conc0[4] = 1.5e-5
  conc0[5] = 1.0                                         # Water: pure-liquid activity ≡ 1 (buffered; enters Q as 1)
  conc0[6] = so4 / 96000.0
  conc0[8] = 1e-7
  conc0[15] = ALK / 1000000.0

  @printf("[%d/%d] SO₄²⁻ = %g mg/L (N=%d strains, deterministic) ... ",
    i, length(SO4_LEVELS), so4, N)
  flush(stdout)
  t0 = time()
  C, T, solver_metadata = full_simulate(
    ps,
    TMAX,
    pop0,
    conc0,
    energy0,
    phi0;
    verbose=false,
    thermodynamic_feedback=true,
    return_metadata=true,
  )
  elapsed = time() - t0
  println(@sprintf("done in %.1f s", elapsed))

  final_pop = C[end, 1:N]
  final_conc = C[end, (N+1):(N+M)]
  final_energy = C[end, (N+M+1):(N+M+N)]
  final_phi = C[end, (N+M+N+1):(N+M+2N)]

  seed = 0
  schema_version = Analysis.OUTPUT_SCHEMA_VERSION
  model_revision = Analysis.MODEL_REVISION
  input_sha256 = Analysis.model_input_hashes()
  thermodynamic_feedback = true
  JLD2.@save joinpath(paras_dir, "anmbr_parameters_run001.jld2") ps seed pop0 conc0 energy0 phi0 schema_version model_revision input_sha256 thermodynamic_feedback solver_metadata
  JLD2.@save joinpath(traj_dir, "anmbr_trajectory_run001.jld2") C T seed final_pop final_conc final_energy final_phi schema_version model_revision input_sha256 thermodynamic_feedback solver_metadata

  open(joinpath(run_dir, "config.toml"), "w") do io
    TOML.print(io, Dict(
      "repeats" => 1,
      "kc" => KC,
      "KS" => KS,
      "kr" => KR,
      "Tmax" => TMAX,
      "initial_population" => INIT_POP,
      "initial_energy" => INIT_ENERGY,
      "initial_fraction" => INIT_FRAC,
      "glucose_mM" => GLUCOSE_MM,
      "alkalinity_mgCaCO3_L" => ALK,
      "sulfate_mg_L" => so4,
      "seed" => 0,
      "seed_sequence" => [0],
      "output_root" => string(ABLATION_ROOT),
      "save_trajectories" => true,
      "thermodynamic_feedback" => true,
      "schema_version" => Analysis.OUTPUT_SCHEMA_VERSION,
      "model_revision" => Analysis.MODEL_REVISION,
      "input_sha256" => Analysis.model_input_hashes(),
      "julia_version" => string(VERSION),
      "generated_at_utc" => string(now(UTC)),
      "adm1_like" => true,
      "guild_sizes" => GUILD_SIZES,
    ))
  end

  CSV.write(joinpath(run_dir, "summary.csv"), DataFrame(
    run=[1], seed=[0],
    survivors=[count(>(1e-9), final_pop)],
    extinct=[count(<=(1e-9), final_pop)],
    final_time=[T[end]], runtime_seconds=[elapsed],
    biomass=[sum(final_pop)],
    mean_energy=[mean(final_energy)],
    mean_fraction=[mean(final_phi)],
    CH4_M=[final_conc[14]],
    SO4_M=[final_conc[6]],
    HS_M=[final_conc[7]],
    solver_retcode=[solver_metadata.retcode],
    solver_success=[solver_metadata.successful],
    reached_tmax=[solver_metadata.reached_tmax],
  ))

  open(joinpath(run_dir, "COMPLETED.toml"), "w") do io
    TOML.print(io, Dict(
      "schema_version" => Analysis.OUTPUT_SCHEMA_VERSION,
      "model_revision" => Analysis.MODEL_REVISION,
      "input_sha256" => Analysis.model_input_hashes(),
      "completed_at_utc" => string(now(UTC)),
      "runs" => 1,
      "all_solver_success" => solver_metadata.successful,
      "all_reached_tmax" => solver_metadata.reached_tmax,
    ))
  end

  push!(created_dirs, run_dir)
end

println()
println("="^72)
println("ADM1-like sweep complete. Created dirs:")
for d in created_dirs
  println("  ", d)
end

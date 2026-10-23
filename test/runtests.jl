using Analysis
using Random
using Test

include(joinpath(@__DIR__, "..", "src", "simulation_registry.jl"))

const METABOLITE_COMPOSITION = [
  (6.0, 12.0, 6.0, 0.0, 0.0),   # glucose
  (1.0, 1.0, 3.0, 0.0, -1.0),   # bicarbonate
  (2.0, 3.0, 2.0, 0.0, -1.0),   # acetate
  (1.0, 0.0, 2.0, 0.0, 0.0),    # carbon dioxide
  (0.0, 2.0, 1.0, 0.0, 0.0),    # water
  (0.0, 0.0, 4.0, 1.0, -2.0),   # sulfate
  (0.0, 1.0, 0.0, 1.0, -1.0),   # bisulfide
  (0.0, 1.0, 0.0, 0.0, 1.0),    # proton
  (3.0, 5.0, 3.0, 0.0, -1.0),   # lactate
  (4.0, 7.0, 2.0, 0.0, -1.0),   # butyrate
  (3.0, 5.0, 2.0, 0.0, -1.0),   # propionate
  (5.0, 9.0, 2.0, 0.0, -1.0),   # valerate
  (0.0, 2.0, 0.0, 0.0, 0.0),    # hydrogen
  (1.0, 4.0, 0.0, 0.0, 0.0),    # methane
  (1.0, 0.0, 3.0, 0.0, -2.0),   # carbonate
]

function reaction_balance(reaction)
  balance = zeros(5)
  for (index, coefficient) in zip(reaction.Rct, reaction.stoichR)
    balance .-= coefficient .* collect(METABOLITE_COMPOSITION[index])
  end
  for (index, coefficient) in zip(reaction.Prd, reaction.stoichP)
    balance .+= coefficient .* collect(METABOLITE_COMPOSITION[index])
  end
  return balance
end

@testset "ANMBR model invariants" begin
  Random.seed!(1)
  ps = Analysis.initialise_anmbr(
    10.0,
    0.001375,
    10.0;
    κ_sulfate=1.875e-8,
    verbose=false,
  )

  @test (ps.N, ps.M, ps.O) == (60, 15, 22)
  for reaction in ps.reacs
    @test reaction_balance(reaction) ≈ zeros(5) atol=1e-8
  end

  pop = fill(1e5, ps.N)
  conc = zeros(ps.M)
  conc[[1, 2, 4, 5, 6, 8, 15]] .= [5e-3, 1.2e-2, 1.5e-5, 1.0, 30 / 96000, 1e-7, 1.2e-3]
  energy = fill(1e5, ps.N)
  phi = fill(0.128, ps.N)
  state = [pop; conc; energy; phi]
  original = copy(state)
  derivative = similar(state)
  rates = zeros(ps.N, ps.O)
  Analysis.AssemblyCook.full_dynamics!(derivative, state, ps, rates, 0.0, true)
  @test state == original
  @test all(isfinite, derivative)

  mic = ps.mics[1]
  reaction = ps.reacs[mic.Reacs[1]]
  substrates = fill(1e-6, length(reaction.Rct))
  products = fill(1e3, length(reaction.Prd))
  enzyme = Analysis.AssemblyCook.Eα(0.128, mic, 1)
  rate_on = Analysis.qs(substrates, products, enzyme, 1, mic, ps.T, reaction)
  rate_off = Analysis.qs(
    substrates,
    products,
    enzyme,
    1,
    mic,
    ps.T,
    reaction;
    thermodynamic_feedback=false,
  )
  @test rate_off >= rate_on >= 0.0

  coupled_14_21 = Analysis.AssemblyCook.compute_coupled_dG0(14, 21, ps.reacs)
  @test coupled_14_21 ≈ ps.reacs[14].ΔG0 + 0.25 * ps.reacs[21].ΔG0
  coupled_15_21 = Analysis.AssemblyCook.compute_coupled_dG0(15, 21, ps.reacs)
  @test coupled_15_21 ≈ ps.reacs[15].ΔG0 + 0.5 * ps.reacs[21].ΔG0
end

@testset "Simulation metadata and solver mode" begin
  config = Analysis.anmbrRunConfig(
    repeats=1,
    output_root="simulations/test-output",
    thermodynamic_feedback=false,
  )
  @test !config.thermodynamic_feedback
  @test isabspath(config.output_root)

  ps = Analysis.initialise_syn_trio(10.0, 0.001375, 10.0)
  C, T, metadata = Analysis.full_simulate(
    ps,
    10.0,
    fill(1e3, ps.N),
    [2.5e-3, 2.5e-3, 0.0],
    fill(1e5, ps.N),
    fill(0.128, ps.N);
    verbose=false,
    thermodynamic_feedback=false,
    return_metadata=true,
  )
  @test metadata.successful
  @test metadata.reached_tmax
  @test !metadata.thermodynamic_feedback
  @test size(C, 1) == length(T)
end

@testset "Run registry metadata" begin
  @test "full_ensemble_2026_08_01" in list_run_sets()
  @test "theta_ablation_2026_08_03" in list_run_sets()
  @test run_set_feedback("full_ensemble_2026_08_01")
  @test !run_set_feedback("theta_ablation_2026_08_03")
  @test get_run_set("theta_ablation_2026_08_03")["control_set"] == "full_ensemble_2026_08_01"
  ensemble = get_run_set("full_ensemble_2026_08_01")
  @test Float64[entry["so4_mg_L"] for entry in ensemble["dirs"]] == [10.0, 30.0, 100.0, 300.0, 1000.0]
  @test run_file_number("anmbr_parameters_run042.jld2") == 42
  @test Analysis.OUTPUT_SCHEMA_VERSION == 3
  @test manual_source_id([joinpath("root-a", "same")]) ==
        manual_source_id([joinpath("root-a", "same")])
  @test manual_source_id([joinpath("root-a", "same")]) !=
        manual_source_id([joinpath("root-b", "same")])
end

@testset "Schema-3 run-directory validation" begin
  mktempdir() do run_dir
    paras_dir = joinpath(run_dir, "paras")
    trajectories_dir = joinpath(run_dir, "trajectories")
    mkpath(paras_dir)
    mkpath(trajectories_dir)

    seed = 7
    schema_version = Analysis.OUTPUT_SCHEMA_VERSION
    model_revision = Analysis.MODEL_REVISION
    input_sha256 = Analysis.model_input_hashes()
    thermodynamic_feedback = true
    solver_metadata = (
      successful=true,
      reached_tmax=true,
      retcode="Success",
    )
    ps = Analysis.initialise_syn_trio(10.0, 0.001375, 10.0)
    pop0 = fill(1e3, ps.N)
    conc0 = zeros(ps.M)
    energy0 = fill(1e5, ps.N)
    phi0 = fill(0.128, ps.N)
    C = reshape([pop0; conc0; energy0; phi0], 1, :)
    T = [0.0]
    final_pop = copy(pop0)
    final_conc = copy(conc0)
    final_energy = copy(energy0)
    final_phi = copy(phi0)

    parameter_path = joinpath(paras_dir, "anmbr_parameters_run001.jld2")
    trajectory_path = joinpath(trajectories_dir, "anmbr_trajectory_run001.jld2")
    JLD2.@save parameter_path ps seed pop0 conc0 energy0 phi0 schema_version model_revision input_sha256 thermodynamic_feedback solver_metadata
    JLD2.@save trajectory_path C T seed final_pop final_conc final_energy final_phi schema_version model_revision input_sha256 thermodynamic_feedback solver_metadata
    CSV.write(
      joinpath(run_dir, "summary.csv"),
      DataFrame(
        run=[1],
        seed=[seed],
        solver_retcode=["Success"],
        solver_success=[true],
        reached_tmax=[true],
      ),
    )
    open(joinpath(run_dir, "config.toml"), "w") do io
      TOML.print(io, Dict(
        "schema_version" => schema_version,
        "model_revision" => model_revision,
        "input_sha256" => input_sha256,
        "thermodynamic_feedback" => thermodynamic_feedback,
        "repeats" => 1,
        "seed_sequence" => [seed],
        "julia_version" => string(VERSION),
        "generated_at_utc" => "2026-07-22T00:00:00Z",
        "sulfate_mg_L" => 30.0,
      ))
    end
    open(joinpath(run_dir, "COMPLETED.toml"), "w") do io
      TOML.print(io, Dict(
        "schema_version" => schema_version,
        "model_revision" => model_revision,
        "input_sha256" => input_sha256,
        "completed_at_utc" => "2026-07-22T00:00:01Z",
        "runs" => 1,
        "all_solver_success" => true,
        "all_reached_tmax" => true,
      ))
    end

    @test validate_run_directory(run_dir)
    @test run_directory_feedback(run_dir)
    @test run_files_for_seed(run_dir, seed).run == 1
    @test validate_paired_community(ps, ps; seed=seed)
    @test validate_paired_run_directories([run_dir, run_dir]; seeds=[seed])
    @test validate_seed_condition_grid(
      DataFrame(seed=[seed], SO4_mgL=[30.0]),
      [run_dir],
    )

    manifest = write_source_manifest(
      run_dir;
      script=@__FILE__,
      manifest_name="test.source.toml",
    )
    @test manifest["script_sha256"] == file_sha256(@__FILE__)
  end
end

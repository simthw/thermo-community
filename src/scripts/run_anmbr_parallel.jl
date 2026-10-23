#!/usr/bin/env julia
"""
Run ANMBR simulation in parallel

Usage:
  julia --project=Analysis -p 8 Analysis/src/scripts/run_anmbr_parallel.jl [options]

Where -p 8 specifies using 8 worker processes (can be adjusted based on CPU cores)

Command line arguments:
  --repeats N        Number of runs (default: 50)
  --tmax TIME        Maximum simulation time (default: 2.0e6)
  --kc VALUE         Catalytic rate constant (default: 10.0)
  --seed N           Random seed (default: 1)
  --output_root DIR  Raw simulation directory root (default: Analysis/simulations)
  --no_save_trajectories  Do not save full trajectories
  --no_thermodynamic_feedback  Use the pure Monod rate law
"""

using Distributed
using ArgParse

# ========== Load module on all worker processes (must be at top level) ==========
@everywhere begin
  using Pkg
  Pkg.activate(joinpath(@__DIR__, "..", ".."))
  using Analysis
end

# ========== Main program ==========

function build_parser()
  settings = ArgParseSettings(
    description="Run ANMBR simulation in parallel"
  )

  @add_arg_table settings begin
    "--repeats"
    help = "Number of independent simulations"
    arg_type = Int
    default = 50
    "--kc"
    help = "Catalytic rate constant"
    arg_type = Float64
    default = 10.0
    "--KS"
    help = "Half-saturation constant"
    arg_type = Float64
    default = (1 / 4) * 5.5e-3
    "--kr"
    help = "Reversibility constant"
    arg_type = Float64
    default = 10.0
    "--tmax"
    help = "Maximum simulation time"
    arg_type = Float64
    default = 2.0e6
    "--initial_population"
    help = "Initial population size"
    arg_type = Float64
    default = 1e5
    "--initial_energy"
    help = "Initial energy reserve"
    arg_type = Float64
    default = 1.0e5
    "--initial_fraction"
    help = "Initial ribosome fraction"
    arg_type = Float64
    default = 0.128
    "--glucose_mM"
    help = "Initial glucose concentration (mM)"
    arg_type = Float64
    default = 5.0
    "--alkalinity"
    help = "Alkalinity (mg CaCO3/L)"
    arg_type = Float64
    default = 1200.0
    "--sulfate"
    help = "Sulfate concentration (mg/L)"
    arg_type = Float64
    default = 30.0
    "--seed"
    help = "Random seed base"
    arg_type = Int
    default = 1
    "--output_root"
    help = "Output directory root path"
    arg_type = String
    default = ""
    "--no_save_trajectories"
    help = "Do not save full trajectories"
    action = :store_true
    "--no_thermodynamic_feedback"
    help = "Disable theta feedback and use the pure Monod rate law"
    action = :store_true
  end

  return settings
end

function main()
  parser = build_parser()
  args = parse_args(parser; as_symbols=false)

  println("=" ^ 60)
  println("ANMBR Parallel Simulation")
  println("=" ^ 60)
  println("Number of workers: ", Distributed.nworkers())
  println("Number of simulations: ", args["repeats"])
  println("=" ^ 60)

  # Build configuration
  output_root = isempty(strip(get(args, "output_root", ""))) ? Analysis.DEFAULT_OUTPUT_ROOT : args["output_root"]

  config = Analysis.anmbrRunConfig(
    repeats=args["repeats"],
    kc=args["kc"],
    KS=args["KS"],
    kr=args["kr"],
    Tmax=args["tmax"],
    initial_population=args["initial_population"],
    initial_energy=args["initial_energy"],
    initial_fraction=args["initial_fraction"],
    glucose_mM=args["glucose_mM"],
    alkalinity_mgCaCO3_L=args["alkalinity"],
    sulfate_mg_L=args["sulfate"],
    seed=args["seed"],
    output_root=output_root,
    save_trajectories=!args["no_save_trajectories"],
    thermodynamic_feedback=!args["no_thermodynamic_feedback"],
  )

  # Run parallel simulation
  run_dir = Analysis.run_anmbr_analysis_parallel(config)

  println("\n" * "=" ^ 60)
  println("Simulation complete!")
  println("Results directory: $run_dir")
  println("=" ^ 60)
end

main()

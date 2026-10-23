#!/usr/bin/env julia
using ArgParse
using Analysis

function build_parser()
  settings = ArgParseSettings(
    description="Reproduce Assembly anmbr simulations using the Analysis module."
  )

  @add_arg_table settings begin
    "--repeats"
    help = "Number of independent simulations to run."
    arg_type = Int
    default = 5
    "--kc"
    help = "Mean catalytic rate constant passed to initialise_anmbr."
    arg_type = Float64
    default = 10.0
    "--KS"
    help = "Mean half-saturation constant passed to initialise_anmbr."
    arg_type = Float64
    default = (1 / 4) * 5.5e-3
    "--kr"
    help = "Mean reversibility constant passed to initialise_anmbr."
    arg_type = Float64
    default = 10.0
    "--tmax"
    help = "Maximum integration time supplied to full_simulate."
    arg_type = Float64
    default = 2.0e6
    "--initial_population"
    help = "Initial cell count for each strain."
    arg_type = Float64
    default = 1e5
    "--initial_energy"
    help = "Initial ATP-equivalent reserves for each strain."
    arg_type = Float64
    default = 1.0e5
    "--initial_fraction"
    help = "Initial ribosome fraction for each strain."
    arg_type = Float64
    default = 0.128
    "--glucose_mM"
    help = "Initial glucose concentration in mM."
    arg_type = Float64
    default = 5.0
    "--alkalinity"
    help = "Alkalinity as mg CaCO3/L (typical: 1000-1500)."
    arg_type = Float64
    default = 1200.0
    "--sulfate"
    help = "Sulfate concentration in mg/L (baseline: 30, high: 960)."
    arg_type = Float64
    default = 30.0
    "--seed"
    help = "Base RNG seed used to seed each run sequentially."
    arg_type = Int
    default = 1
    "--output_root"
    help = "Directory root where outputs will be stored."
    arg_type = String
    "--no_save_trajectories"
    help = "Disable saving full state trajectories (only summary/parameters saved)."
    action = :store_true
    "--no_thermodynamic_feedback"
    help = "Disable theta feedback and use the pure Monod rate law."
    action = :store_true
  end

  return settings
end

function derive_config(args::Dict{String,Any})
  output_root_arg = get(args, "output_root", nothing)
  output_root = isnothing(output_root_arg) || isempty(strip(String(output_root_arg))) ? Analysis.DEFAULT_OUTPUT_ROOT : String(output_root_arg)
  save_traj = !get(args, "no_save_trajectories", false)

  return Analysis.anmbrRunConfig(
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
    save_trajectories=save_traj,
    thermodynamic_feedback=!get(args, "no_thermodynamic_feedback", false),
  )
end

function main()
  parser = build_parser()
  args = parse_args(parser; as_symbols=false)
  println("Parsed arguments:")
  for (key, value) in args
    println("  $key => $value")
  end
  config = derive_config(args)
  run_dir = Analysis.run_anmbr_analysis(config)
  println("Results stored in: ", run_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end

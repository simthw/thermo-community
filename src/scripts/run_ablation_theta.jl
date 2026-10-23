#!/usr/bin/env julia
# θ ablation — parallel SO4 sweep with thermodynamic feedback DISABLED.

# Mechanism: sets `thermodynamic_feedback=false` in the serialized run config,
# making the rate law collapse to pure Monod (q = kc·E·S / (KS + S)).

# Parallelization: spawns N worker processes (default 20), loads Analysis on each,
# then dispatches per-seed simulations through `run_anmbr_analysis_parallel`.

# Outputs: `Analysis/simulations/ablation_theta/{timestamp}/` per SO4 level.

# Usage (defaults: 100 seeds × 5 SO4, 20 workers):
#   julia --project=Analysis Analysis/src/scripts/run_ablation_theta.jl

# Env vars (PowerShell):
#   $env:ABLATION_THETA_SEEDS = "10"       # pilot run with 10 seeds
#   $env:ABLATION_THETA_WORKERS = "8"      # use 8 cores instead of 20
#   $env:ABLATION_THETA_SO4 = "10,1000"    # custom SO4 list

using Distributed
using Dates
using Printf

# ===========================================================================
# 1. Read config from env vars BEFORE spawning workers
# ===========================================================================
const SEEDS = parse(Int, get(ENV, "ABLATION_THETA_SEEDS", "100"))
const N_WORKERS = parse(Int, get(ENV, "ABLATION_THETA_WORKERS", "20"))
const SO4_LEVELS = parse.(Float64, split(get(ENV, "ABLATION_THETA_SO4", "10,30,100,300,1000"), ","))

# ===========================================================================
# 2. Spawn workers — they start cold and need Analysis loaded individually
# ===========================================================================
@printf("Spawning %d worker processes ...\n", N_WORKERS)
addprocs(N_WORKERS)
@printf("Workers ready: master + %d workers (total processes: %d)\n",
  nworkers(), nprocs())
println()

# ===========================================================================
# 3. Load Analysis on master + all workers
# ===========================================================================
println("Loading Analysis on all processes (precompilation may take ~1 min the first time) ...")
@everywhere using Analysis
println("Analysis loaded everywhere.")
println()

# ===========================================================================
# 4. Confirm the explicit, serializable model mode
# ===========================================================================
println("θ ablation mode: thermodynamic_feedback=false (pure Monod)")
println()

# ===========================================================================
# 5. Sweep configuration
# ===========================================================================
const ABLATION_ROOT = joinpath(Analysis.DEFAULT_OUTPUT_ROOT, "ablation_theta")
isdir(ABLATION_ROOT) || mkpath(ABLATION_ROOT)

@printf("Seeds per SO₄²⁻: %d\n", SEEDS)
@printf("SO₄²⁻ levels (mg/L): %s\n", string(SO4_LEVELS))
@printf("Output root: %s\n", ABLATION_ROOT)
@printf("Parallel workers: %d\n", N_WORKERS)
println()

# ===========================================================================
# 6. Run sweep — call the parallel entry point
# ===========================================================================
created_dirs = String[]
for (i, so4) in enumerate(SO4_LEVELS)
  println("="^72)
  @printf("[%d/%d]  SO₄²⁻ = %g mg/L\n", i, length(SO4_LEVELS), so4)
  println("="^72)

  cfg = anmbrRunConfig(
    repeats=SEEDS,
    sulfate_mg_L=so4,
    output_root=ABLATION_ROOT,
    thermodynamic_feedback=false,
  )
  out_dir = run_anmbr_analysis_parallel(cfg)
  push!(created_dirs, out_dir)
end

println()
println("="^72)
println("θ ablation sweep complete.")
println("="^72)
println("Output dirs (ordered by SO₄²⁻ ascending):")
for d in created_dirs
  println("  ", d)
end

module AssemblyCook

using DifferentialEquations
using Distributions
using LinearAlgebra
using SciMLBase
using Statistics
using SymPy
using JLD2
import ..Analysis: Rgas, ΔGATP, NA, Keq

include("fullParameters.jl")
include("make_full.jl")
include("make_full_robust.jl")
include("simulate_full.jl")
include("thermodynamics.jl")
include("analytic.jl")
include("interactions.jl")
include("robustness.jl")

# Export from fullParameters.jl
export Reaction, make_Reaction
export MicrobeP, make_MicrobeP
export FullParameters, make_FullParameters, extinction

# Export from make_full.jl
export initialise_anmbr, initialise, initialise_syn_trio, initialise_pol_duo, initialise_adm1_like

# Export from make_full_robust.jl
export options_titles, initialise_np, initialise_Pb, initialise_ϕH, initialise_sat

# Export from simulate_full.jl
export full_simulate, full_simulate_syn
export θ, qs

# Export from thermodynamics.jl
export calculate_reaction_gibbs, analyze_thermodynamic_feasibility
export calculate_actual_ΔG, reaction_entropy_production, total_entropy_production
export get_reaction_energy_sum

# Export from analytic.jl
export Force, nForce, Jacobian, nJacobian

# Export from interactions.jl
export quantify_ints, analyse_ints, InteractionResult

# Export from robustness.jl
export robustness_test, batch_robustness_test, RobustnessResult

end # module

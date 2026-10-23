# Script that tests robustness of communities by generating communities
# under different parameterisations
# Note: JLD2 and Distributed are imported in AssemblyCook.jl

# Export core functions
export robustness_test, RobustnessResult

"""
    RobustnessResult

Struct to hold robustness test results.
"""
struct RobustnessResult
    option::Int64  # Option number (0-6)
    option_name::String  # Option name
    N::Int64  # Number of strains
    survivors::Int64  # Number of surviving strains
    extinct::Int64  # Number of extinct strains
    final_time::Float64  # Final simulation time
    runtime_seconds::Float64  # Runtime in seconds
    ps::FullParameters  # Final parameter set (with extinct removed)
end

"""
    robustness_test(N::Int64, M::Int64, O::Int64, Rl::Int64, Ru::Int64,
    kc::Float64, KS::Float64, kr::Float64, syn::Bool, μrange::Float64,
    opt::Int64, Tmax::Float64=1e8; nworkers::Int64=1)

Test robustness under different parameterisations.

# Arguments
- `N::Int64`: Number of strains
- `M::Int64`: Number of metabolites
- `O::Int64`: Number of reactions
- `Rl::Int64`: Lower bound on reactions per strain
- `Ru::Int64`: Upper bound on reactions per strain
- `kc::Float64`: Catalytic rate constant
- `KS::Float64`: Half-saturation constant
- `kr::Float64`: Reversibility constant
- `syn::Bool`: Whether syntrophy is enabled
- `μrange::Float64`: Free energy range
- `opt::Int64`: Robustness option (0-6)
- `Tmax::Float64`: Maximum simulation time (default: 1e8)
- `nworkers::Int64`: Number of parallel workers (default: 1)

# Returns
- `RobustnessResult`: Test results
"""
function robustness_test(N::Int64, M::Int64, O::Int64, Rl::Int64, Ru::Int64,
    kc::Float64, KS::Float64, kr::Float64, syn::Bool, μrange::Float64,
    opt::Int64, Tmax::Float64=1e8; nworkers::Int64=1)
    
    # Check inputs
    @assert Rl >= 1 "Lower bound must be at least 1"
    @assert Ru >= Rl "Upper bound cannot be lower than lower bound"
    @assert N >= 1 "Number of strains must be positive"
    @assert opt >= 0 && opt <= 6 "Option must be between 0 and 6"
    
    # Get option name
    option_name = options_titles(opt)
    
    # Initial conditions
    pop = ones(N)
    conc = zeros(M)
    as = 1e5 * ones(N)
    ϕs = 0.128 * ones(N)
    
    # Generate parameter set based on option
    if opt == 0
        # Standard case
        ps = initialise(N, M, Rl, Ru, kc, KS, kr, syn; κ_sulfate=1.875e-8)
    elseif opt == 1
        # Option 1 changes average mass of metabolic protein
        np = 600
        ps = initialise_np(N, M, O, Rl, Ru, kc, KS, kr, syn, μrange, np)
    elseif opt == 2
        # Option 2 changes average bound fraction of ribosomes
        Pb = 0.9
        ps = initialise_Pb(N, M, O, Rl, Ru, kc, KS, kr, syn, μrange, Pb)
    elseif opt == 3
        # Option 3 changes size of housekeeping fraction
        ϕH = 0.7
        ps = initialise_ϕH(N, M, O, Rl, Ru, kc, KS, kr, syn, μrange, ϕH)
    elseif opt == 4
        # Option 4 makes Kγ bigger than KΩ
        Kγ = 1e9
        KΩ = 5e8
        ps = initialise_sat(N, M, O, Rl, Ru, kc, KS, kr, syn, μrange, Kγ, KΩ)
    elseif opt == 5
        # Option 5 makes Kγ and KΩ both bigger
        Kγ = 5e9
        KΩ = 1e10
        ps = initialise_sat(N, M, O, Rl, Ru, kc, KS, kr, syn, μrange, Kγ, KΩ)
    else
        # Option 6 makes Kγ and KΩ both smaller
        Kγ = 5e7
        KΩ = 1e8
        ps = initialise_sat(N, M, O, Rl, Ru, kc, KS, kr, syn, μrange, Kγ, KΩ)
    end
    
    # Run simulation
    t_start = time()
    C, T = full_simulate(ps, Tmax, pop, conc, as, ϕs; verbose=false)
    elapsed = time() - t_start
    
    # Establish which microbes are extinct (below e-10 threshold)
    ext = (C[end, 1:N] .<= 1e-10)
    survivors = count(x -> !x, ext)
    extinct = count(x -> x, ext)
    
    # Remove extinct species from parameter set
    ps_final = extinction(ps, ext)
    
    return RobustnessResult(
        opt,
        option_name,
        N,
        survivors,
        extinct,
        T[end],
        elapsed,
        ps_final
    )
end

"""
    batch_robustness_test(N::Int64, M::Int64, O::Int64, Rl::Int64, Ru::Int64,
    kc::Float64, KS::Float64, kr::Float64, syn::Bool, μrange::Float64,
    options::Vector{Int64}=collect(0:6); Tmax::Float64=1e8)

Run robustness tests for multiple options in batch.

# Arguments
- `N::Int64`: Number of strains
- `M::Int64`: Number of metabolites
- `O::Int64`: Number of reactions
- `Rl::Int64`: Lower bound on reactions per strain
- `Ru::Int64`: Upper bound on reactions per strain
- `kc::Float64`: Catalytic rate constant
- `KS::Float64`: Half-saturation constant
- `kr::Float64`: Reversibility constant
- `syn::Bool`: Whether syntrophy is enabled
- `μrange::Float64`: Free energy range
- `options::Vector{Int64}`: List of options to test (default: 0-6)
- `Tmax::Float64`: Maximum simulation time (default: 1e8)

# Returns
- `Vector{RobustnessResult}`: Results for all options
"""
function batch_robustness_test(N::Int64, M::Int64, O::Int64, Rl::Int64, Ru::Int64,
    kc::Float64, KS::Float64, kr::Float64, syn::Bool, μrange::Float64,
    options::Vector{Int64}=collect(0:6); Tmax::Float64=1e8)
    
    results = RobustnessResult[]
    
    for opt in options
        println("Testing option $opt: $(options_titles(opt))")
        result = robustness_test(N, M, O, Rl, Ru, kc, KS, kr, syn, μrange, opt, Tmax)
        push!(results, result)
        println("  Survivors: $(result.survivors)/$(result.N)")
        println("  Runtime: $(round(result.runtime_seconds, digits=2))s")
    end
    
    return results
end

# Script that makes the parameters needed for simulation of the full proteome model
export initialise, initialise_anmbr, initialise_syn_trio, initialise_pol_duo, initialise_adm1_like

using CSV
using DataFrames

# --- Thermodynamics Calculation Functions ---

"""
    _load_metabolite_energies()

Load metabolite Gibbs free energy of formation data from CSV file.

Returns a tuple of (ΔGf0_dict, ΔHf0_dict) where:
- ΔGf0_dict: Dictionary mapping chemical formula to ΔGf0 (J/mol)
- ΔHf0_dict: Dictionary mapping chemical formula to ΔHf0 (J/mol) or nothing if unavailable
"""
function _load_metabolite_energies()
  csv_path = joinpath(@__DIR__, "..", "data", "metabolite_energies.csv")

  if !isfile(csv_path)
    error("Metabolite energies file not found: $csv_path")
  end

  df = CSV.read(csv_path, DataFrame)

  ΔGf0_dict = Dict{String,Float64}()
  ΔHf0_dict = Dict{String,Union{Float64,Nothing}}()

  for row in eachrow(df)
    formula = strip(row.formula)
    ΔGf0_dict[formula] = row.dGf0_kJ_per_mol * 1000.0  # Convert kJ/mol to J/mol

    # Handle enthalpy data (may be missing or "unknown")
    if !ismissing(row.dHf0_kJ_per_mol) && !(row.dHf0_kJ_per_mol isa AbstractString)
      ΔHf0_dict[formula] = row.dHf0_kJ_per_mol * 1000.0
    else
      ΔHf0_dict[formula] = nothing
    end
  end

  return ΔGf0_dict, ΔHf0_dict
end

"""
    _calculate_reaction_ΔG0(rct_formulas, rct_coeffs, prd_formulas, prd_coeffs, ΔGf0_dict)

Calculate the standard Gibbs free energy change for a reaction.

# Arguments
- `rct_formulas::Vector{String}`: Reactant chemical formulas
- `rct_coeffs::Vector{Float64}`: Reactant stoichiometric coefficients
- `prd_formulas::Vector{String}`: Product chemical formulas
- `prd_coeffs::Vector{Float64}`: Product stoichiometric coefficients
- `ΔGf0_dict::Dict{String, Float64}`: Dictionary of standard Gibbs free energies of formation

# Returns
- `Float64`: ΔG0 in J/mol

# Formula
ΔG0 = Σ(ν_P × ΔGf0_P) - Σ(ν_R × ΔGf0_R)
"""
function _calculate_reaction_ΔG0(
  rct_formulas::Vector{String},
  rct_coeffs::Vector{Float64},
  prd_formulas::Vector{String},
  prd_coeffs::Vector{Float64},
  ΔGf0_dict::Dict{String,Float64}
)
  # Calculate products contribution
  ΔG0_products = 0.0
  for i in 1:length(prd_formulas)
    formula = prd_formulas[i]
    if !haskey(ΔGf0_dict, formula)
      error("Unknown metabolite in products: $formula")
    end
    ΔG0_products += prd_coeffs[i] * ΔGf0_dict[formula]
  end

  # Calculate reactants contribution
  ΔG0_reactants = 0.0
  for i in 1:length(rct_formulas)
    formula = rct_formulas[i]
    if !haskey(ΔGf0_dict, formula)
      error("Unknown metabolite in reactants: $formula")
    end
    ΔG0_reactants += rct_coeffs[i] * ΔGf0_dict[formula]
  end

  return ΔG0_products - ΔG0_reactants
end

# --- Public Constants ---
const DEFAULT_T = 310.0
const DEFAULT_MC = 10^8
const DEFAULT_γm = 1260.0 / 60.0
const DEFAULT_n = [7459, 300, 300]
const DEFAULT_d = 6.0e-5  # Matches Assembly project
const DEFAULT_ρ = 29.0
const DEFAULT_Kγ = 5e8
const DEFAULT_Pb = 0.7
const DEFAULT_ϕH = 0.45
const DEFAULT_KΩ = 1e9
const DEFAULT_fd = log(100) / log(2)
const DEFAULT_κ_OFFSET = 3.3e-7  # Metabolite supply rate, matches Assembly project

# Bicarbonate (HCO3-) supply rate to maintain alkalinity buffer
# Target alkalinity: 1000-1500 mg CaCO3/L ≈ 0.01-0.015 M HCO3-
# At steady state: C_ss = κ/δ, so κ = δ * C_target
# With δ = 6e-5 and C_target = 0.012 M (1200 mg CaCO3/L):
# κ_HCO3 = 6e-5 * 0.012 = 7.2e-7 M/s
const DEFAULT_κ_HCO3 = 7.2e-7  # Bicarbonate supply rate for alkalinity buffering

# this function randomly assigns a unique set of reaction identifiers to a microbe, selecting
# a number of reactions between a specified minimum and maximum from a total pool, ensuring no duplicates.
function choose_reactions(O::Int64, Rl::Int64, Ru::Int64)
  @assert Ru <= O "Strain cannot have more reactions than the number of possible reactions"
  R = rand(Rl:Ru)
  Reacs = zeros(Int64, R)
  Reacs[1] = rand(1:O)
  for i in 2:R
    good = false
    while good == false
      r = rand(1:O)
      if r ∉ Reacs[1:(i-1)]
        Reacs[i] = r
        good = true
      end
    end
  end
  return (R, Reacs)
end

# function that assigns reactions based on strain type
function choose_reactions_from_type(available::Vector{Int64}, Rl::Int64, Ru::Int64)
  O = length(available)
  @assert Ru <= O "Strain cannot have more reactions than available for its type"
  R = rand(Rl:Ru)
  Reacs = zeros(Int64, R)
  for i in 1:R
    good = false
    while !good
      idx = rand(1:O)
      r = available[idx]
      if r ∉ Reacs[1:(i-1)]
        Reacs[i] = r
        good = true
      end
    end
  end
  return (R, Reacs)
end

# Hydrogen production and consumption reaction IDs
const H2_PRODUCTION_REACTIONS = [14, 15, 16, 17, 18, 19, 20] # VFAs → H2

# H2 consumption reactions split by substrate availability:

# RELIABLE: reaction 21 (hydrogenotrophic methanogenesis, 4H2 + CO2 → CH4)
#   - substrate CO2/HCO3- is always available (alkalinity buffer maintained by κ[2])
#   - valid basis for static η assignment at initialisation time

# CONDITIONAL: reaction 10 (sulfate reduction, H2 + SO42- → HS-)
#   - substrate SO42- supply is intentionally low (κ[6] = δ·[SO4]/96000;
#     ≈1.9e-8 M/s at the 30 mg/L baseline)
#   - sulfate may be transiently depleted at runtime
#   - its effect on H2 partial pressure is already captured dynamically by θ = Q/Keq:
#     when SO42- → 0, θ(rxn 10) → 1, sulfate reducers stall, H2 accumulates,
#     which in turn raises θ for H2-producing reactions and slows them down
#   - must NOT be used for η assignment; doing so would inflate Keq for H2 producers
#     even when sulfate is absent, giving them an unrealistic thermodynamic advantage
const H2_CONSUMPTION_REACTIONS_RELIABLE = [21]     # always-available partners: used for η assignment
const H2_CONSUMPTION_REACTIONS_CONDITIONAL = [10]    # substrate-limited partners: handled by θ at runtime
const H2_CONSUMPTION_REACTIONS = vcat(H2_CONSUMPTION_REACTIONS_RELIABLE,
  H2_CONSUMPTION_REACTIONS_CONDITIONAL)

"""
 compute_coupled_dG0(h2_prod_id, h2_cons_id, reacs)

Compute combined ΔG0 for coupled H2 production-consumption system.
This represents the overall thermodynamics when syntrophic partners consume H2.
"""
function compute_coupled_dG0(h2_prod_id::Int, h2_cons_id::Int, reacs::Array{Reaction,1})
  producer = reacs[h2_prod_id]
  consumer = reacs[h2_cons_id]
  producer_h2 = findfirst(==(13), producer.Prd)
  consumer_h2 = findfirst(==(13), consumer.Rct)
  producer_h2 === nothing && error("Reaction $h2_prod_id does not produce H2")
  consumer_h2 === nothing && error("Reaction $h2_cons_id does not consume H2")

  # Scale the consumer reaction so its H2 demand equals the producer's H2
  # yield. Adding unscaled reactions is only valid when both coefficients match.
  coupling_ratio = producer.stoichP[producer_h2] / consumer.stoichR[consumer_h2]
  dG_combined = producer.ΔG0 + coupling_ratio * consumer.ΔG0

  return dG_combined
end

"""
 find_best_syntrophic_partner(h2_prod_id, available_reacs, reacs)

Find the best *reliable* H2-consuming reaction for static η assignment.

Only H2_CONSUMPTION_REACTIONS_RELIABLE (reaction 21, hydrogenotrophic methanogenesis)
is considered. Conditional partners such as sulfate reduction (reaction 11) are excluded
because their substrate (SO42-) may be depleted at runtime; their thermodynamic effect on
H2 partial pressure is already captured dynamically by θ = Q/Keq during simulation.

Returns (best_partner_id, combined_dG0), or (nothing, Inf) if no reliable partner exists.
"""
function find_best_syntrophic_partner(h2_prod_id::Int, available_reacs::Array{Int64,1},
  reacs::Array{Reaction,1})
  best_partner = nothing
  best_dG = Inf

  for cons_id in H2_CONSUMPTION_REACTIONS_RELIABLE   # reliable partners only
    if cons_id ∈ available_reacs
      dG_combined = compute_coupled_dG0(h2_prod_id, cons_id, reacs)
      if dG_combined < best_dG
        best_dG = dG_combined
        best_partner = cons_id
      end
    end
  end

  return best_partner, best_dG
end

# function to choose uniformly mixed values for η
function choose_η_mix(reacs::Array{Reaction,1}, Reacs::Array{Int64,1},
  T::Float64, syn::Bool; verbose::Bool=false)
  η = zeros(length(Reacs))
  ηl = 1 / 3 # minimum ATP yield for exergonic reactions
  mratio = syn == true ? 1e-2 : 1e5

  # Check if the reaction pool contains any reliable H2-consuming reactions (e.g. reaction 21,
  # hydrogenotrophic methanogenesis). This is a community-level property: H2 producers can
  # benefit from syntrophy even though the H2 consumer lives in a different microbe.
  #
  # Note: `syn` is NOT used here. In the original Assembly model, syn merely controls mratio
  # (the minimum equilibrium product/substrate ratio), which sets the upper bound on η for
  # exergonic reactions. It has no meaning related to community composition.
  has_h2_consumer = any(r ∈ H2_CONSUMPTION_REACTIONS_RELIABLE for r in 1:length(reacs))

  for i in eachindex(η)
    I = Reacs[i]
    dG = reacs[I].ΔG0

    if dG > 0.0
      # Endergonic reaction (e.g., H₂ production from VFAs)

      if I ∈ H2_PRODUCTION_REACTIONS && has_h2_consumer
        # Syntrophic H2 production: search all 22 reactions for the best H2-consuming partner,
        # not just reactions carried by this microbe. H2 is an inter-species intermediate —
        # the consumer lives in a different cell, so the partner search must span the full
        # community reaction pool.
        all_reaction_ids = collect(1:length(reacs))
        best_partner, dG_combined = find_best_syntrophic_partner(I, all_reaction_ids, reacs)

        if best_partner !== nothing && dG_combined < 0.0
          # Coupled system is exergonic! Allow H2 production with reduced η
          # The organism gets partial benefit from the coupled system
          # η_eff reflects the "syntrophic bonus"

          # Calculate effective η based on combined free energy
          # The more negative dG_combined, the higher the potential η
          ηh_coupled = -(dG_combined + Rgas * T * log(mratio)) / ΔGATP
          ηh_coupled = max(ηh_coupled, ηl)

          # Use lower bound for syntrophic reactions (more conservative)
          # Range: [ηl, ηh_coupled *0.5] - syntrophic organisms get less direct benefit
          η_max = ηh_coupled * 0.5
          η[i] = ηl + (η_max - ηl) * rand() * 0.5 # Conservative random value
          η[i] = max(η[i], ηl)

          if verbose
            println(" Syntrophic H₂ production (reaction $I): coupled with reaction $best_partner, " *
                    "combined ΔG0 = $(round(dG_combined / 1000, digits=2)) kJ/mol, η = $(round(η[i], digits=3))")
          end
        else
          # No favorable partner found, treat as standard endergonic
          η[i] = 0.0
        end
      else
        # Standard endergonic reaction (non-H2 producing or no H2 consumers)
        # Set η =0 so that Keq = exp(-ΔG0/RT) reflects pure thermodynamic equilibrium;
        # the organism gains energy from the actual ΔG at in-situ concentrations,
        # which is computed in full_dynamics! via the ΔG_actual formula.
        η[i] = 0.0
      end
    else
      # Exergonic reaction: standard treatment
      ηh = -(dG + Rgas * T * log(mratio)) / (ΔGATP)
      ηh = max(ηh, ηl)
      η[i] = (ηh - ηl) * rand() + ηl
      η[i] = max(η[i], ηl)
    end
  end
  return η
end

# function to generate fix set of reaction for simple test modules. The steps between metabolites are fixed.
function fix_reactions(O::Int64, M::Int64, μrange::Float64, T::Float64)
  @assert O == 2 * M - 3 "Miscalculated the number of reactions expected"
  RP = zeros(Int64, O, 2)
  ΔG = zeros(O)
  dG = -μrange / (M - 1)
  for i in 1:O
    if i % 2 != 0
      RP[i, 1] = ceil(i / 2)
      RP[i, 2] = ceil(i / 2) + 1
      ΔG[i] = dG
    else
      RP[i, 1] = ceil(i / 2)
      RP[i, 2] = ceil(i / 2) + 2
      ΔG[i] = 2 * dG
    end
  end
  return (RP, ΔG)
end

# function to take in average kinetic parameters and return randomised vectors of them
function kin_rand(kc::Float64, KS::Float64, kr::Float64, R::Int64)
  d = Normal()
  rs = rand(d, R, 3)
  kcs = kc * (2.0 .^ rs[:, 1])
  KSs = KS * (2.0 .^ rs[:, 2])
  krs = kr * (2.0 .^ rs[:, 3])
  return (kcs, KSs, krs)
end

"""
    _make_microbe_params_with_κ(M::Int64, κ_values::Vector{Float64})

Create a NamedTuple with default microbe parameters and custom κ values.
"""
function _make_microbe_params_with_κ(M::Int64, κ_values::Vector{Float64})
  return (
    T=DEFAULT_T,
    MC=DEFAULT_MC,
    γm=DEFAULT_γm,
    n=copy(DEFAULT_n),
    d=DEFAULT_d,
    ρ=DEFAULT_ρ,
    Kγ=DEFAULT_Kγ,
    Pb=DEFAULT_Pb,
    ϕH=DEFAULT_ϕH,
    KΩ=DEFAULT_KΩ,
    fd=DEFAULT_fd,
    δ=DEFAULT_d * ones(M),
    κ=copy(κ_values),
  )
end

"""
    _build_microbe(params, Reacs, η, kcs, KSs, krs)

Helper function to construct a MicrobeP object with default parameters.
"""
function _build_microbe(
  params::NamedTuple,
  Reacs::Vector{Int64},
  η::Vector{Float64},
  kcs::Vector{Float64},
  KSs::Vector{Float64},
  krs::Vector{Float64}
)
  R = length(Reacs)
  ϕP = rand(R)
  # Normalize so sum(ϕP) = 1.0, ensuring all available protein is allocated to metabolic enzymes
  # In Assembly project, make_MicrobeP requires sum(ϕP) ≈ 1.0
  ϕP = ϕP ./ sum(ϕP)
  return make_MicrobeP(
    params.MC, params.γm, params.ρ, params.Kγ, params.Pb, params.d,
    params.ϕH, params.KΩ, params.fd, R, Reacs, η, kcs, KSs, krs, params.n, ϕP
  )
end

# --- Reaction Parsing Helpers ---
const METABOLITE_DICT = Dict(
  "C6H12O6(aq)" => 1, "HCO3-(aq)" => 2, "CH3COO-(aq)" => 3, "CO2(g)" => 4, "H2O(l)" => 5,
  "SO42-(aq)" => 6, "HS-(aq)" => 7, "H+(aq)" => 8, "CH3CH(OH)COO-(aq)" => 9,
  "CH3CH2CH2COO-(aq)" => 10, "CH3CH2COO-(aq)" => 11, "CH3CH2CH2CH2COO-(aq)" => 12,
  "H2(g)" => 13, "CH4(g)" => 14, "CO32-(aq)" => 15
)

# Note: ΔG0_VALUES dictionary has been removed as ΔG0 is now calculated
# automatically from metabolite Gibbs free energies in metabolite_energies.csv

# Reaction type mapping (based on reactions.dat definition)
# 1-6:   Acidogenic fermentation (Glucose → Acetate/Lactate/Butyrate/Valerate)
# 7-11:  Sulfate reduction (VFAs + SO42- → HS-)
# 12-13: Acetogenesis (Homoacetogenesis, not methanogenesis!)
# 14-20: Hydrogen production (VFAs → H2)
# 21-22: Methanogenesis (H2/Acetate → CH4) ← Only these two produce CH4
const REACTION_TYPE_DICT = Dict(
  1 => "acidifier", 2 => "acidifier", 3 => "acidifier", 4 => "acidifier",
  5 => "acidifier", 6 => "acidifier", 7 => "sulfate_reducer", 8 => "sulfate_reducer",
  9 => "sulfate_reducer", 10 => "sulfate_reducer", 11 => "sulfate_reducer",
  12 => "acetogen", 13 => "acetogen", 14 => "hydrogen_producer",
  15 => "hydrogen_producer", 16 => "hydrogen_producer", 17 => "hydrogen_producer",
  18 => "hydrogen_producer", 19 => "hydrogen_producer", 20 => "hydrogen_producer",
  21 => "methanogen", 22 => "methanogen"
)

function _parse_term(term::AbstractString)
  term = strip(term)
  if occursin(r"^\d", term)
    parts = split(term, " ", limit=2)
    coeff = parse(Float64, parts[1])
    formula = strip(parts[2])
  else
    coeff = 1.0
    formula = term
  end
  # Keep the state marker (aq), (g), (l) as it's part of the metabolite key
  return coeff, formula
end

function _parse_reactions_dat()
  # Load metabolite Gibbs free energy data
  ΔGf0_dict, ΔHf0_dict = _load_metabolite_energies()

  lines = readlines(joinpath(@__DIR__, "..", "data", "reactions.dat"))
  reacs = Array{Reaction,1}()
  reaction_id = 1

  for line in lines
    if startswith(strip(line), "#") || isempty(strip(line))
      continue
    end

    parts = split(line, ": ")
    if length(parts) < 2
      continue
    end

    equation = strip(parts[2])

    # Remove everything after # (comment) - but we no longer parse ΔG0 from comments
    if occursin("#", equation)
      comment_idx = findfirst('#', equation)
      equation = strip(equation[1:(comment_idx-1)])
    end

    # Split reactants and products (using -> instead of → for better compatibility)
    if occursin("->", equation)
      reactants_str, products_str = split(equation, "->")
    elseif occursin("→", equation)
      reactants_str, products_str = split(equation, "→")
    else
      error("Invalid reaction format: missing '->' or '→' in equation")
    end

    reactants_str = strip(reactants_str)
    products_str = strip(products_str)

    # Parse reactants
    rct_indices = Int64[]
    rct_coeffs = Float64[]
    rct_formulas = String[]

    for term in split(reactants_str, " + ")
      coeff, formula = _parse_term(String(term))
      if haskey(METABOLITE_DICT, formula)
        push!(rct_indices, METABOLITE_DICT[formula])
        push!(rct_coeffs, coeff)
        push!(rct_formulas, formula)
      else
        error("Unknown metabolite in reactants: $formula")
      end
    end

    # Parse products
    prd_indices = Int64[]
    prd_coeffs = Float64[]
    prd_formulas = String[]

    for term in split(products_str, " + ")
      coeff, formula = _parse_term(String(term))
      if haskey(METABOLITE_DICT, formula)
        push!(prd_indices, METABOLITE_DICT[formula])
        push!(prd_coeffs, coeff)
        push!(prd_formulas, formula)
      else
        error("Unknown metabolite in products: $formula")
      end
    end

    # Calculate ΔG0 from metabolite Gibbs free energies
    ΔG0 = _calculate_reaction_ΔG0(rct_formulas, rct_coeffs,
      prd_formulas, prd_coeffs,
      ΔGf0_dict)

    # Create Reaction object
    push!(reacs, Reaction(reaction_id, rct_indices, prd_indices,
      rct_coeffs, prd_coeffs, ΔG0))
    reaction_id += 1
  end

  return reacs, length(reacs)
end

function _build_type_to_reactions_map(type_dict::Dict{Int,String}, O::Int64)
  type_to_reactions = Dict{String,Vector{Int}}()
  # Anaerobic digestion food chain:
  # acidifier(fermentation) → hydrogen_producer(H2 production) → methanogen(methanogenesis)
  #                                           ↘ sulfate_reducer(sulfate reduction, competes with methanogen for H2)
  # acetogen: homoacetogenic bacteria, consumes H2 to produce acetate
  strain_types = ["acidifier", "hydrogen_producer", "methanogen", "sulfate_reducer", "acetogen"]
  for (id, typ) in type_dict
    if id > O
      continue
    end
    if !haskey(type_to_reactions, typ)
      type_to_reactions[typ] = Int[]
    end
    push!(type_to_reactions[typ], id)
  end
  return type_to_reactions, strain_types
end

# Baseline strain type proportions for mesophilic AnMBR with low-sulfate influent.
# Order: [acidifier, hydrogen_producer, methanogen, sulfate_reducer, acetogen]
const STRAIN_TYPE_RATIOS = [0.35, 0.20, 0.25, 0.10, 0.10]  # sum should be 1.0

function _get_strain_type_idx(i::Int, N::Int)
  cumsum_ratios = cumsum(STRAIN_TYPE_RATIOS)
  normalized_pos = (i - 0.5) / N
  for (idx, threshold) in enumerate(cumsum_ratios)
    normalized_pos <= threshold && return idx
  end
  return length(STRAIN_TYPE_RATIOS)
end

# function to generate parameter set for the model with inhibition
function initialise(N::Int64, M::Int64, Rl::Int64, Ru::Int64, kc::Float64,
  KS::Float64, kr::Float64, syn::Bool; κ_sulfate::Float64,
  κ_HCO3::Float64=DEFAULT_κ_HCO3, verbose::Bool=false)
  κ = zeros(M)
  κ[1] = DEFAULT_κ_OFFSET     # Glucose supply (C6H12O6) - main carbon source feed
  κ[2] = κ_HCO3               # Bicarbonate supply (HCO3-); scales with configured alkalinity so the buffer pool matches the target
  κ[6] = κ_sulfate             # Sulfate (SO42-) supply rate, configurable for different influent concentrations
  # Note: Acetate, H2 and other intermediate products are not supplied externally, produced by microbial metabolism
  params = _make_microbe_params_with_κ(M, κ)

  reacs, O = _parse_reactions_dat()  # Dynamically parse reaction pool
  type_to_reactions, strain_types = _build_type_to_reactions_map(REACTION_TYPE_DICT, O)

  mics = Array{MicrobeP,1}(undef, N)
  for i in 1:N
    type_idx = _get_strain_type_idx(i, N)
    type = strain_types[type_idx]
    available = type_to_reactions[type]
    effective_Ru = min(Ru, length(available))
    R, Reacs = choose_reactions_from_type(available, Rl, effective_Ru)
    kcs, KSs, krs = kin_rand(kc, KS, kr, R)
    η = choose_η_mix(reacs, Reacs, params.T, syn; verbose=verbose)
    mics[i] = _build_microbe(params, Reacs, η, kcs, KSs, krs)
  end

  O = length(reacs)
  ps = make_FullParameters(N, M, O, params.T, params.κ, params.δ, reacs, mics)
  return (ps)
end

"""
    initialise_anmbr(kc, KS, kr; κ_sulfate, verbose=false)

Generate a parameter set for simulating a 60-strain AnMBR community with 15 metabolites.
Strain types are allocated according to fixed baseline proportions in `STRAIN_TYPE_RATIOS`.

# Arguments
- `kc::Float64` : Variability parameter for catalytic constants
- `KS::Float64` : Variability parameter for half-saturation constants
- `kr::Float64` : Variability parameter for inhibition constants
- `κ_sulfate::Float64` : Sulfate supply rate (M/s).
    Compute from influent concentration: κ_sulfate = δ × (sulfate_mg_L / 96000)
- `κ_HCO3::Float64` : Bicarbonate supply rate (M/s). Controls the HCO3- buffer pool;
    pH is held at 7 under the assumption that this buffer is adequate.
    Compute from configured alkalinity: κ_HCO3 = δ × (alkalinity_mgCaCO3_L / 100000).

# Returns
- `FullParameters` : Parameter set for AnMBR simulation
"""
function initialise_anmbr(kc::Float64, KS::Float64, kr::Float64;
  κ_sulfate::Float64, κ_HCO3::Float64=DEFAULT_κ_HCO3, verbose::Bool=false)
  N = 60
  M = 15
  Rl = 1
  Ru = 5
  syn = true
  return initialise(N, M, Rl, Ru, kc, KS, kr, syn;
    κ_sulfate=κ_sulfate, κ_HCO3=κ_HCO3, verbose=verbose)
end

"""
    compute_deterministic_η(reacs, Reacs, T, syn; verbose=false)

Deterministic counterpart to `choose_η_mix`. Within each branch of `choose_η_mix`,
the random `rand()` call is replaced by the midpoint of the permissible η range.
Endergonic reactions without a viable syntrophic partner still receive η=0.
"""
function compute_deterministic_η(reacs::Array{Reaction,1}, Reacs::Array{Int64,1},
  T::Float64, syn::Bool; verbose::Bool=false)
  η = zeros(length(Reacs))
  ηl = 1 / 3
  mratio = syn == true ? 1e-2 : 1e5
  has_h2_consumer = any(r ∈ H2_CONSUMPTION_REACTIONS_RELIABLE for r in 1:length(reacs))

  for i in eachindex(η)
    I = Reacs[i]
    dG = reacs[I].ΔG0
    if dG > 0.0
      if I ∈ H2_PRODUCTION_REACTIONS && has_h2_consumer
        all_reaction_ids = collect(1:length(reacs))
        best_partner, dG_combined = find_best_syntrophic_partner(I, all_reaction_ids, reacs)
        if best_partner !== nothing && dG_combined < 0.0
          ηh_coupled = -(dG_combined + Rgas * T * log(mratio)) / ΔGATP
          ηh_coupled = max(ηh_coupled, ηl)
          η_max = ηh_coupled * 0.5
          # original draws from [ηl, ηl + (η_max - ηl)*0.5] via `rand()*0.5`;
          # deterministic midpoint:
          η[i] = ηl + (η_max - ηl) * 0.25
          η[i] = max(η[i], ηl)
          if verbose
            println("  Syntrophic H₂ producer rx $I (det.): coupled with rx $best_partner, " *
                    "combined ΔG0 = $(round(dG_combined/1000, digits=2)) kJ/mol, η = $(round(η[i], digits=3))")
          end
        else
          η[i] = 0.0
        end
      else
        η[i] = 0.0
      end
    else
      ηh = -(dG + Rgas * T * log(mratio)) / ΔGATP
      ηh = max(ηh, ηl)
      η[i] = (ηh + ηl) / 2    # midpoint of [ηl, ηh]
      η[i] = max(η[i], ηl)
    end
  end
  return η
end

"""
    initialise_adm1_like(kc, KS, kr; κ_sulfate, κ_HCO3=DEFAULT_κ_HCO3, verbose=false)

Deterministic five-guild baseline: ONE strain per guild, each carrying the
FULL set of guild reactions, with deterministic (median/midpoint) parameters.

Construction:
- 5 strains, one per `strain_types` entry from `_build_type_to_reactions_map`:
  acidifier, hydrogen_producer, methanogen, sulfate_reducer, acetogen
- Each strain's Reacs = all reactions in its guild (no random subset; `sort`-ed
  for deterministic ordering)
- `kcs = fill(kc, R)` — no log-normal multiplier
- `KSs = fill(KS, R)`, `krs = fill(kr, R)` — same
- `η` from `compute_deterministic_η` (midpoint of `choose_η_mix`'s range)
- `ϕP = fill(1/R, R)` — uniform proteome allocation (no random Dirichlet)

The same `M=15` metabolite system and thermodynamic feedback theta are kept.
This bundled contrast changes strain count, reaction coverage, kinetic
variability, eta, and proteome allocation together; it is not a single-factor
diversity intervention and is not an implementation of ADM1 equations.

Use this as a deterministic reference curve, not as an orthogonal causal proof
that within-guild diversity alone is necessary for response multiplicity.

# Arguments
- `kc::Float64`, `KS::Float64`, `kr::Float64`: deterministic kinetic constants
- `κ_sulfate::Float64`: sulfate supply rate (M/s)
- `κ_HCO3::Float64`: bicarbonate supply rate (M/s)

# Returns
- `FullParameters` with `N=5`
"""
function initialise_adm1_like(kc::Float64, KS::Float64, kr::Float64;
  κ_sulfate::Float64, κ_HCO3::Float64=DEFAULT_κ_HCO3, verbose::Bool=false)
  M = 15
  κ = zeros(M)
  κ[1] = DEFAULT_κ_OFFSET
  κ[2] = κ_HCO3
  κ[6] = κ_sulfate
  params = _make_microbe_params_with_κ(M, κ)

  reacs, O = _parse_reactions_dat()
  type_to_reactions, strain_types = _build_type_to_reactions_map(REACTION_TYPE_DICT, O)

  N = length(strain_types)
  mics = Array{MicrobeP,1}(undef, N)
  syn = true

  for i in 1:N
    type = strain_types[i]
    Reacs = sort(copy(type_to_reactions[type]))
    R = length(Reacs)

    kcs = fill(kc, R)
    KSs = fill(KS, R)
    krs = fill(kr, R)
    η = compute_deterministic_η(reacs, Reacs, params.T, syn; verbose=verbose)
    ϕP = fill(1.0 / R, R)

    mics[i] = make_MicrobeP(
      params.MC, params.γm, params.ρ, params.Kγ, params.Pb, params.d,
      params.ϕH, params.KΩ, params.fd, R, Reacs, η, kcs, KSs, krs, params.n, ϕP
    )
    if verbose
      println("  ADM1-like strain $i ($type): R=$R, Reacs=$Reacs, η=$(round.(η, digits=3))")
    end
  end

  return make_FullParameters(N, M, O, params.T, params.κ, params.δ, reacs, mics)
end

# function to generate set of three strains that interact thermodynamically
function initialise_syn_trio(kc::Float64, KS::Float64, kr::Float64)
  N = 3
  M = 3
  O = 3
  R = 1

  κ = zeros(M)
  κ[1] = DEFAULT_κ_OFFSET
  params = _make_microbe_params_with_κ(M, κ)

  μrange = 5e6 * (M / 25)
  RP, ΔG = fix_reactions(O, M, μrange, params.T)

  reacs = Array{Reaction,1}(undef, O)
  for i in 1:O
    reacs[i] = Reaction(i, [RP[i, 1]], [RP[i, 2]], [1.0], [1.0], ΔG[i])
  end

  mics = Array{MicrobeP,1}(undef, N)
  for i in 1:N
    if i == 1
      Reacs, η = [1], [3.0]
    elseif i == 2
      Reacs, η = [1], [3.85]
    elseif i == 3
      Reacs, η = [3], [3.5]
    else
      error("too many strains provided, need a trio")
    end
    kcs, KSs, krs = [kc], [KS], [kr]
    mics[i] = _build_microbe(params, Reacs, η, kcs, KSs, krs)
  end

  ps = make_FullParameters(N, M, O, params.T, params.κ, params.δ, reacs, mics)
  return (ps)
end

# function to generate set of two strains that interact thermodynamically
function initialise_pol_duo(kc::Float64, KS::Float64, kr::Float64)
  N = 2
  M = 3
  O = 3
  R = 1

  κ = zeros(M)
  κ[1] = DEFAULT_κ_OFFSET
  κ[2] = DEFAULT_κ_OFFSET
  params = _make_microbe_params_with_κ(M, κ)

  μrange = 5e6 * (M / 25)
  RP, ΔG = fix_reactions(O, M, μrange, params.T)

  reacs = Array{Reaction,1}(undef, O)
  for i in 1:O
    reacs[i] = Reaction(i, [RP[i, 1]], [RP[i, 2]], [1.0], [1.0], ΔG[i])
  end

  mics = Array{MicrobeP,1}(undef, N)
  for i in 1:N
    if i == 1
      Reacs, η = [3], [4.00]
    elseif i == 2
      Reacs, η = [2], [6.00]
    else
      error("too many strains provided, need a duo")
    end
    kcs, KSs, krs = [kc], [KS], [kr]
    mics[i] = _build_microbe(params, Reacs, η, kcs, KSs, krs)
  end

  ps = make_FullParameters(N, M, O, params.T, params.κ, params.δ, reacs, mics)
  return (ps)
end

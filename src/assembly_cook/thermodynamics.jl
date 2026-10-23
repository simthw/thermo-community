export calculate_reaction_gibbs, analyze_thermodynamic_feasibility,
  calculate_actual_ΔG, reaction_entropy_production, total_entropy_production,
  get_reaction_energy_sum, dissipation

"""
    calculate_reaction_gibbs(rxn::Reaction, conc::Vector{Float64}, T::Float64)

Calculate the actual Gibbs free energy change (ΔG) for a reaction given metabolite concentrations.
Returns ΔG in kJ/mol.
ΔG = ΔG0 + RT * ln(Q)
"""
function calculate_reaction_gibbs(rxn::Reaction, conc::Vector{Float64}, T::Float64)
  # Calculate Reaction Quotient Q
  # Q = (Π [P]^ν_P) / (Π [S]^ν_S)

  numerator = 1.0
  for (i, prod_idx) in enumerate(rxn.Prd)
    # Ensure concentration is non-zero to avoid log(0)
    c = max(conc[prod_idx], 1e-12)
    numerator *= c^rxn.stoichP[i]
  end

  denominator = 1.0
  for (i, rct_idx) in enumerate(rxn.Rct)
    c = max(conc[rct_idx], 1e-12)
    denominator *= c^rxn.stoichR[i]
  end

  Q = numerator / denominator

  # ΔG0 in Reaction struct is in J/mol
  ΔG_J = rxn.ΔG0 + Rgas * T * log(Q)

  return ΔG_J / 1000.0 # Return in kJ/mol for readability
end

"""
    calculate_reaction_gibbs(rxn::Reaction, conc_dict::Dict{Int, Float64}, T::Float64, M::Int64)
Overload for dictionary input.
"""
function calculate_reaction_gibbs(rxn::Reaction, conc_dict::Dict{Int,Float64}, T::Float64, M::Int64)
  conc_vec = zeros(M)
  for (k, v) in conc_dict
    if k <= M
      conc_vec[k] = v
    end
  end
  return calculate_reaction_gibbs(rxn, conc_vec, T)
end

"""
    analyze_thermodynamic_feasibility(ps::FullParameters, conc::Vector{Float64})
Print a report of thermodynamic feasibility for all reactions.
"""
function analyze_thermodynamic_feasibility(ps::FullParameters, conc::Vector{Float64})
  println("\n--- Thermodynamic Feasibility Analysis (T=$(ps.T) K) ---")
  println(rpad("Rxn ID", 8) * rpad("Type", 10) * rpad("ΔG (kJ/mol)", 15) * "Status")

  for rxn in ps.reacs
    dG = calculate_reaction_gibbs(rxn, conc, ps.T)
    status = dG < 0 ? "Feasible" : "Infeasible (Blocked)"

    # Determine simplified type based on reactants/products
    rtype = "Unknown"
    if 13 in rxn.Prd
      rtype = "H₂ Prod"
    end
    if 14 in rxn.Prd
      rtype = "Methanogen"
    end
    if 7 in rxn.Prd
      rtype = "Sulf Red"
    end

    println(rpad(string(rxn.ID), 8) * rpad(rtype, 10) * rpad(string(round(dG, digits=2)), 15) * status)
  end
  println("------------------------------------------------------\n")
end

"""
    calculate_actual_ΔG(S, P, T, ΔG0, stoichR, stoichP)

Calculate actual Gibbs free energy change for a reaction.

# Arguments
- `S::Vector{Float64}`: Substrate concentrations
- `P::Vector{Float64}`: Product concentrations
- `T::Float64`: Temperature (K)
- `ΔG0::Float64`: Standard free energy change (J/mol)
- `stoichR::Vector{Float64}`: Substrate stoichiometric coefficients
- `stoichP::Vector{Float64}`: Product stoichiometric coefficients

# Returns
- `Float64`: Actual ΔG (J/mol)
"""
function calculate_actual_ΔG(S::Vector{Float64}, P::Vector{Float64}, T::Float64,
  ΔG0::Float64, stoichR::Vector{Float64}, stoichP::Vector{Float64})
  R = Rgas  # J/(mol·K)

  # Calculate reaction quotient Q
  # Q = (Π [P]^ν_P) / (Π [S]^ν_S)
  Q_numerator = prod(max.(P, 1e-12) .^ stoichP)
  Q_denominator = prod(max.(S, 1e-12) .^ stoichR)
  Q = Q_numerator / Q_denominator

  # ΔG = ΔG0 + RT ln(Q)
  ΔG = ΔG0 + R * T * log(Q)
  return ΔG
end

"""
    reaction_entropy_production(rate, ΔG, T)

Calculate entropy production rate for a single reaction.

# Arguments
- `rate::Float64`: Reaction rate (mol/s)
- `ΔG::Float64`: Actual Gibbs free energy change (J/mol)
- `T::Float64`: Temperature (K)

# Returns
- `Float64`: Entropy production rate (J/(K·s))
"""
function reaction_entropy_production(rate::Float64, ΔG::Float64, T::Float64)
  # The project reports the magnitude of the Clausius term per reaction.
  σ = abs(rate * ΔG) / T
  return σ
end

"""
    total_entropy_production(C::Matrix{Float64}, ps::FullParameters, t_idx::Int=1)

Calculate total ecosystem entropy production rate at a given time step.

This follows the Clausius entropy-production measure used in the current
Analysis project:
    σ = |J · ΔG| / T
    Total σ = Σ |q · ΔG · pop / NA| / T

Where:
- ΔG = ΔG0 + R·T·ln(Q), excluding ATP-capture correction
- q is the per-cell reaction-event rate
- pop is population size (cells)
- NA is Avogadro's number

# Arguments
- `C::Matrix{Float64}`: Full state matrix, rows=time steps, columns=variables
  - C[:, 1:N] = populations
  - C[:, N+1:N+M] = metabolite concentrations
  - C[:, N+M+1:N+M+N] = energy (ATP)
  - C[:, N+M+N+1:N+M+2N] = ribosome fractions
- `ps::FullParameters`: System parameters
- `t_idx::Int`: Time step index (default=1)

# Returns
- `Float64`: Total entropy production rate (J/(K·s))
"""
function total_entropy_production(
  C::Matrix{Float64},
  ps::FullParameters,
  t_idx::Int=1;
  thermodynamic_feedback::Bool=true,
)
  σ_total = 0.0
  T = ps.T
  N = ps.N
  M = ps.M

  # Extract state from specified time step
  C_row = C[t_idx, :]
  pop = C_row[1:N]
  conc = C_row[(N+1):(N+M)]
  ϕR_vec = C_row[(N+M+N+1):(N+M+2N)]

  for i in 1:N  # Loop over each microbe
    # Skip extinct strains
    if pop[i] <= 1e-10
      continue
    end

    mic = ps.mics[i]

    for j in 1:mic.R  # Loop over each reaction
      reac_idx = mic.Reacs[j]
      r = ps.reacs[reac_idx]

      # Get substrate and product concentrations
      S_vec = [conc[rct] for rct in r.Rct]
      P_vec = [conc[prd] for prd in r.Prd]

      ΔG = calculate_actual_ΔG(S_vec, P_vec, T, r.ΔG0, r.stoichR, r.stoichP)

      # Calculate enzyme amount Eα = MC * (1 - ϕR - ϕH) * ϕP[j] / n[2]
      ϕR = ϕR_vec[i]
      E = mic.MC * max(1.0 - ϕR - mic.ϕH, 0.0) * mic.ϕP[j] / mic.n[2]

      q = qs(
        S_vec,
        P_vec,
        E,
        j,
        mic,
        T,
        r;
        thermodynamic_feedback=thermodynamic_feedback,
      )

      if q <= 0.0
        continue
      end

      # Calculate entropy production rate σ = (q · ΔG) / T
      # Note: q is molecular-level rate, convert to molar units
      σ = (q * ΔG * pop[i] / NA) / T
      σ_total += abs(σ)
    end
  end

  return σ_total
end

"""
    dissipation(C_row, N, M, ps)

Calculate free energy dissipation rate for the ecosystem at a single time step.

This function follows the Assembly project's thermodynamic formulation:
    Fd = -(ΔG0 + R·T·ln(Q) + η·ΔGATP)
    Total dissipation = Σ q · Fd · pop / NA

Where Fd is the free energy dissipated per reaction event (total ΔG minus ATP captured).
Positive Fd means net dissipation (energy lost as heat).

# Arguments
- `C_row::Vector{Float64}`: State vector at one time step (length = N + M + N + N)
  - C_row[1:N] = populations
  - C_row[N+1:N+M] = metabolite concentrations
  - C_row[N+M+1:N+M+N] = energy (ATP)
  - C_row[N+M+N+1:N+M+2N] = ribosome fractions
- `N::Int`: Number of strains
- `M::Int`: Number of metabolites
- `ps::FullParameters`: System parameters

# Returns
- `dsp_total::Float64`: Total ecosystem dissipation rate (J/s)
- `dsp_by_strain::Vector{Float64}`: Dissipation rate per strain (J/s)
"""
function dissipation(C_row::Vector{Float64}, N::Int, M::Int, ps::FullParameters)
  Temp = ps.T
  pop = C_row[1:N]
  conc = C_row[(N+1):(N+M)]
  ϕR_vec = C_row[(N+M+N+1):(N+M+2N)]

  dsp_total = 0.0
  dsp_by_strain = zeros(N)

  for i in 1:N
    # Skip extinct strains (population effectively zero)
    if pop[i] <= 0.0
      continue
    end

    mic = ps.mics[i]
    for j in 1:mic.R
      reac_idx = mic.Reacs[j]
      rxn = ps.reacs[reac_idx]

      # Get substrate and product concentrations
      S_vec = [conc[rct] for rct in rxn.Rct]
      P_vec = [conc[prd] for prd in rxn.Prd]

      # Skip if any product is zero (log(0) undefined)
      if any(P_vec .<= 0.0)
        continue
      end

      # Reaction quotient Q = Π P^stoichP / Π S^stoichR
      S_safe = max.(S_vec, 1e-15)
      Q_val = prod(P_vec .^ rxn.stoichP) / prod(S_safe .^ rxn.stoichR)

      # Free energy dissipated per reaction event (J/mol):
      # Fd = -(ΔG0 + RT·ln(Q) + η·ΔGATP)
      # = -(actual ΔG) - η·ΔGATP
      # Positive Fd means net dissipation (energy lost as heat)
      Fd = -(rxn.ΔG0 + Rgas * Temp * log(Q_val) + mic.η[j] * ΔGATP)

      # Enzyme amount Eα = MC * (1 - ϕR - ϕH) * ϕP[j] / n[2]
      ϕR = ϕR_vec[i]
      E = mic.MC * max(1.0 - ϕR - mic.ϕH, 0.0) * mic.ϕP[j] / mic.n[2]

      # Thermodynamic inhibition θ
      Keq_val = exp((-rxn.ΔG0 - mic.η[j] * ΔGATP) / (Rgas * Temp))
      θ_val = Q_val / Keq_val

      # Reaction rate per cell (Monod kinetics with thermodynamic inhibition)
      S_lim = minimum(max.(S_vec, 0.0) ./ rxn.stoichR)
      if S_lim <= 0.0
        continue
      end
      q = mic.kc[j] * E * S_lim * max(1.0 - θ_val, 0.0) /
          (mic.KS[j] + S_lim * (1.0 + mic.kr[j] * θ_val))
      q = max(q, 0.0)

      if q <= 0.0
        continue
      end

      # Accumulate dissipation: q * Fd * population (molecule-level → /NA for moles)
      d = q * Fd * pop[i] / NA
      dsp_total += d
      dsp_by_strain[i] += d
    end
  end

  return dsp_total, dsp_by_strain
end

"""
    get_reaction_energy_sum(ps::FullParameters)

Calculate sum of standard Gibbs free energies for all reactions in the system.

# Arguments
- `ps::FullParameters`: System parameters

# Returns
- `Float64`: Sum of all ΔG0 values (kJ/mol)
"""
function get_reaction_energy_sum(ps::FullParameters)
  total_ΔG0 = 0.0
  for rxn in ps.reacs
    total_ΔG0 += rxn.ΔG0 / 1000.0  # Convert to kJ/mol
  end
  return total_ΔG0
end

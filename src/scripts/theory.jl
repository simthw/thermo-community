#!/usr/bin/env julia
#region theory.jl - Non-equilibrium thermodynamic steady-state analysis
# Using actual ANMBR reaction data for analysis
# Based on reactions.dat and metabolite_energies.csv
# Reference: Analysis method from 2step.jl
#endregion

# Add project path to load path
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Analysis
include(joinpath(@__DIR__, "..", "figures", "publication_theme.jl"))

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FIGURE_ROOT = joinpath(PROJECT_ROOT, "results", "paper", "theory")

"""
    load_reactions_from_dat()

Load reaction data from reactions.dat, return Reaction vector
"""
function load_reactions_from_dat()
  # Use Analysis module's initialise_anmbr function to get reaction data
  ps = Analysis.initialise_anmbr(10.0, 5.5e-3/4, 10.0; κ_sulfate=1.875e-8)
  return ps.reacs
end

"""
    analyze_reaction_thermodynamics(reac::Reaction)

Analyze thermodynamics of a single reaction
"""
function analyze_reaction_thermodynamics(reac, T::Float64=310.0)
  # Metabolite index to name mapping
  met_names = ["Glucose", "HCO₃⁻", "Acetate", "CO₂", "H₂O",
    "SO₄²⁻", "HS-", "H+", "Lactate", "Butyrate",
    "Propionate", "Valerate", "H₂", "CH₄", "CO₃²⁻"]

  # Get substrate and product names
  substrates = ["$(reac.stoichR[i]) $(met_names[reac.Rct[i]])" for i in 1:length(reac.Rct)]
  products = ["$(reac.stoichP[i]) $(met_names[reac.Prd[i]])" for i in 1:length(reac.Prd)]

  # Calculate equilibrium constant
  η = 1.0  # Assume ATP yield = 1
  Keq_val = Analysis.Keq(T, η, reac.ΔG0)

  return (
    id=reac.ID,
    substrates=substrates,
    products=products,
    ΔG0=reac.ΔG0,
    ΔG0_kJ=reac.ΔG0 / 1000.0,
    Keq=Keq_val
  )
end

"""
    twostep_anaerobic()

Two-step anaerobic digestion model analysis:
- Step 1: Acidifiers ferment glucose to acetate (G → A)
- Step 2: Methanogens convert acetate to methane (A → M)

This model demonstrates:
1. Condition of zero enzyme concentration change at steady state
2. How thermodynamic constraints (θ term) affect reaction flux
3. Energy dissipation distribution between the two steps
"""
function twostep_anaerobic()
  println("=== Two-Step Anaerobic Digestion Model ===")
  println("Step 1: Glucose(G) → Acetate(A) (Acidifier)")
  println("Step 2: Acetate(A) → Methane(M) (Methanogen)")
  println()

  # Select reasonable parameter values
  T = 310.0  # Temperature 310K (37°C, typical for anaerobic digestion)
  cE1 = 1e5  # Acidifier enzyme concentration
  cE2 = 1e5  # Methanogen enzyme concentration
  cG = 0.005 # Glucose concentration 5 mM
  cM = 0.001 # Methane concentration (low)

  # Reaction rate constants
  rk1 = 1000.0  # Acidifier forward
  rk2 = 100.0   # Acidifier reverse
  rk3 = 1000.0  # Methanogen forward
  rk4 = 100.0   # Methanogen reverse

  # Overall equilibrium constant (glucose → methane is strongly exergonic)
  Keq_total = 1e12  # Glucose → Methane equilibrium constant (very large, effectively irreversible)

  # Explore effect of different Keq1 values on steady state
  Keq1_range = [0.01, 0.1, 1.0, 10.0, 100.0]
  Keq2_range = Keq_total ./ Keq1_range

  println("Parameter scan analysis:")
  println("Keq1\t\tKeq2\t\tAcetate(M)\tFlux\tEnergy dissipation(kJ/s)")
  println("-" ^ 75)

  results = []

  for (i, keq1) in enumerate(Keq1_range)
    keq2 = Keq2_range[i]

    # Calculate reverse rate constants
    rK1 = sqrt(rk1 * rk2 / keq1)
    rK2 = rk2 / rK1
    rK3 = sqrt(rk3 * rk4 / keq2)
    rK4 = rk4 / rK3

    # Simplified calculation: assume steady-state acetate concentration
    # Use quasi-steady-state approximation
    A_ss = cG * rk1 * rK3 / (rk4 * rK2 + rk3 * cM)

    # Step 2 steady-state enzyme distribution
    denom2 = A_ss * rk3 + cM * rK4 + rk4 + rK3
    E2_free = cE2 * (rk4 + rK3) / denom2
    E2X = cE2 - E2_free
    flux = E2X * rk4 - E2_free * cM * rK4

    # Calculate energy dissipation (use -ΔG for exothermic)
    Q1 = A_ss / cG  # Step 1 reaction quotient
    Q2 = cM / A_ss  # Step 2 reaction quotient
    ΔG1 = Rgas * T * (log(Q1) - log(keq1))
    ΔG2 = Rgas * T * (log(Q2) - log(keq2))

    # Energy dissipation: energy from exothermic reaction
    energy_dissipated = -((ΔG1 + ΔG2) * flux) / 1000  # Convert to kJ/s

    push!(results, (keq1, keq2, A_ss, flux, energy_dissipated))

    println("$(rpad(round(keq1, digits=2), 10))\t$(rpad(round(keq2, digits=0), 10))\t$(rpad(round(A_ss, digits=6), 12))\t$(rpad(round(flux, digits=2), 12))\t$(round(energy_dissipated, digits=2))")
  end

  println()
  println("Analysis complete!")

  # Visualize results
  try
    plt = publication_figure()
    ax = Axis(plt[1, 1];
      xlabel="Keq1 (step 1 equilibrium constant)",
      ylabel="Steady-state acetate (M)")
    scatterlines!(ax, Keq1_range, [r[3] for r in results];
      marker=:circle, label="Acetate conc.")
    axislegend(ax)
    mkpath(FIGURE_ROOT)
    outpath = joinpath(FIGURE_ROOT, "theory_steady_state.png")
    savefig(outpath, plt)
    println("Figure saved to: $outpath")
  catch e
    println("Warning: Cannot generate figure ($e)")
  end

  return results
end

"""
    thermodynamic_analysis()

Thermodynamic constraint analysis: Calculate θ values for different substrate/product concentrations
"""
function thermodynamic_analysis()
  println("\n=== Thermodynamic Constraint Analysis ===")
  println("θ = Q / Keq")
  println("θ = 0: Reaction proceeds forward, no thermodynamic inhibition")
  println("θ = 1: Reaction at equilibrium, maximum thermodynamic inhibition")
  println("θ > 1: Reaction proceeds in reverse")
  println()

  T = 310.0  # K

  # Use actual reaction ΔG0 (glucose → acetate, typical value ~ -200 kJ/mol)
  ΔG0 = -90000.0  # J/mol (Standard free energy change, including η*ΔGATP correction)
  η = 1.0  # ATP yield

  # Simplified Keq calculation (η·ΔGATP correction omitted)
  Keq_val = exp(-ΔG0 / (Rgas * T))

  println("T = $(T) K")
  println("ΔG0 = $(ΔG0) J/mol")
  println("Keq = $(round(Keq_val, digits=2))")
  println()

  # Calculate θ for different Q values
  Q_values = [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]

  println("Q\t\tθ\t\tReaction status")
  println("-" ^ 50)
  for Q in Q_values
    θ_val = Q / Keq_val
    if θ_val < 0.01
      status = "Forward"
    elseif θ_val < 0.5
      status = "Mild inhibition"
    elseif θ_val < 0.9
      status = "Moderate inhibition"
    elseif θ_val < 1.1
      status = "Near equilibrium"
    elseif θ_val < 10
      status = "Reverse"
    else
      status = "Fully reversed"
    end
    println("$(Q)\t\t$(round(θ_val, digits=6))\t$status")
  end

  # Add actual case calculation
  println()
  println("--- Actual Case ---")
  println("Glucose → Acetate typical ΔG0 ≈ -200 kJ/mol")
  ΔG0_actual = -200000.0  # J/mol
  Keq_actual = exp(-ΔG0_actual / (Rgas * T))
  println("Keq = $(round(Keq_actual, digits=2))")

  # Assume actual concentrations
  G_actual = 0.005  # 5 mM
  A_actual = 0.001  # 1 mM
  Q_actual = A_actual / G_actual
  θ_actual = Q_actual / Keq_actual
  println("Q = A/G = $(Q_actual)")
  println("θ = $(round(θ_actual, digits=6))")
  if θ_actual < 1.0
    println("Reaction can proceed spontaneously (θ < 1)")
  else
    println("Reaction is thermodynamically inhibited (θ ≥ 1)")
  end

  return nothing
end

"""
    entropy_production()

Entropy production analysis: Calculate total entropy production rate of the system
σ = Σ(rate_i * (-ΔG_i)) / T
Negative ΔG indicates exothermic reaction, energy dissipated as entropy
"""
function entropy_production()
  println("\n=== Entropy Production Analysis ===")
  println("σ = Σ(rate_i * (-ΔG_i)) / T")
  println("Negative ΔG indicates exothermic reaction, energy dissipated as entropy to environment")
  println("Non-equilibrium steady state: System continuously produces entropy, energy dissipation required to maintain steady state")
  println()

  # Load actual reaction data from reactions.dat
  reacs = load_reactions_from_dat()

  # Metabolite index to name mapping
  met_names = ["Glucose", "HCO₃⁻", "Acetate", "CO₂", "H₂O",
    "SO₄²⁻", "HS-", "H+", "Lactate", "Butyrate",
    "Propionate", "Valerate", "H₂", "CH₄", "CO₃²⁻"]

  T = 310.0  # K

  # Assumed reaction rate
  rate = 1e-6  # mol/s

  total_entropy = 0.0

  # Group by reaction type
  println("=== Entropy Production Analysis for Major Reactions ===")
  println("Reaction\t\t\t\tΔG (kJ/mol)\tEntropy (J/K·s)")
  println("-" ^ 70)

  # Select several key reactions for analysis
  key_ids = [1, 9, 14, 21, 22]  # Glucose acidification, Butyrate H2 production, Propionate H2 production, Hydrogenotrophic methanogenesis, Acetoclastic methanogenesis

  for rid in key_ids
    reac = reacs[rid]

    # Format reaction
    substrates = join(["$(reac.stoichR[i])$(met_names[reac.Rct[i]])" for i in 1:length(reac.Rct)], " + ")
    products = join(["$(reac.stoichP[i])$(met_names[reac.Prd[i]])" for i in 1:length(reac.Prd)], " + ")

    # Calculate entropy production (use -ΔG because exothermic ΔG < 0)
    dG = reac.ΔG0
    entropy_prod = -rate * dG / T  # dG is negative, so entropy production is positive

    total_entropy += entropy_prod

    # Truncate overly long reaction descriptions
    reaction_str = "$substrates → $products"
    if length(reaction_str) > 25
      reaction_str = first(reaction_str, 22) * "..."
    end

    println("$(rpad(reaction_str, 30))\t$(rpad(round(dG/1000, digits=1), 12))\t$(round(entropy_prod, digits=6))")
  end

  println("-" ^ 70)
  println("Total entropy production rate: $(round(total_entropy, digits=6)) J/(K·s)")
  println()

  # Add explanation
  println("Explanation:")
  println("- Positive entropy production indicates system releases energy to environment")
  println("- Exothermic reactions (ΔG < 0) generate energy, dissipated as entropy")
  println("- At non-equilibrium steady state, system continuously obtains energy from substrates")
  println("- Energy dissipated as entropy, maintaining ordered state (cell growth)")
  println()

  return total_entropy
end

"""
    main()

Main function: Run all theoretical analyses
"""
function main()
  println("Non-equilibrium Thermodynamic Steady-State Analysis")
  println("=" ^ 50)
  println()

  # 1. Two-step anaerobic digestion model
  twostep_anaerobic()

  # 2. Thermodynamic constraint analysis
  thermodynamic_analysis()

  # 3. Entropy production analysis
  entropy_production()

  println("\n" * "=" ^ 50)
  println("Theoretical analysis complete!")
  println()
  println("Key Conclusions:")
  println("1. Steady-state condition: dx/dt = 0 (time derivative of all state variables is zero)")
  println("2. Thermodynamic constraint: θ = Q/Keq controls reaction direction")
  println("3. Non-equilibrium steady state: Continuous mass-energy flow, but constant state")
  println("4. Entropy production: Non-equilibrium systems continuously produce entropy, energy dissipation required to maintain steady state")
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end

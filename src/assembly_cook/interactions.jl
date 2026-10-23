# Script to find and analyse interaction strength and type of strains in simulation
# Note: SymPy, LinearAlgebra, Statistics, JLD2 are imported in AssemblyCook.jl

# Export core functions
export quantify_ints, analyse_ints, InteractionResult

"""
    InteractionResult

Struct to hold interaction analysis results.
"""
struct InteractionResult
  ints::Array{Int64,3}  # Interaction types [N, N, M]
  in_str::Array{Float64,3}  # Interaction strengths [N, N, M]
  net_in::Array{Float64,2}  # Net interactions [N, N]
  Fatp::Array{Float64,2}  # ATP forces [N, M]
  frc::Float64  # Perturbation fraction used
end

"""
    quantify_ints(ps::FullParameters, inf_out::Array{Float64, 1}; frc::Float64=1e-3)

Quantify interaction strength and type for all strains.

# Arguments
- `ps::FullParameters`: System parameters
- `inf_out::Array{Float64, 1}`: Final state vector [populations, metabolites, energies, fractions]
- `frc::Float64`: Fraction to perturb metabolite concentrations (default: 1e-3)

# Returns
- `InteractionResult`: Struct containing interaction analysis results
"""
function quantify_ints(ps::FullParameters, inf_out::Array{Float64,1}; frc::Float64=1e-3)
  # For each strain record metabolites it interacts with
  mts = Array{Array{Int64,1},1}(undef, ps.N)
  for j in 1:(ps.N)
    sbs = reduce(vcat, [ps.reacs[r].Rct for r in ps.mics[j].Reacs])
    pds = reduce(vcat, [ps.reacs[r].Prd for r in ps.mics[j].Reacs])
    # Combine into 1 vector (remove repeats)
    met = unique(vcat(sbs, pds))
    mts[j] = met
  end

  # Preallocate vector of forces (symbolic)
  F = Array{Sym,1}(undef, 3 * ps.N + ps.M)
  F = Force(ps, F)

  # Find initial unperturbed force
  Fi = nForce(F, inf_out, ps)

  # Preallocate array to store interaction types
  ints = zeros(Int64, ps.N, ps.N, ps.M)
  # Preallocate array to store interaction strengths
  in_str = zeros(ps.N, ps.N, ps.M)
  # Preallocate array for forces
  Fatp = zeros(ps.N, ps.M)
  # Preallocate array for net interactions
  net_in = zeros(ps.N, ps.N)

  # Preallocate perturbation vector
  pout = zeros(3 * ps.N + ps.M)
  # Populations, energies, and fractions remain fixed
  pout[1:(ps.N)] = inf_out[1:(ps.N)]
  pout[(ps.N+ps.M+1):(2*ps.N+ps.M)] = inf_out[(ps.N+ps.M+1):(2*ps.N+ps.M)]
  pout[(2*ps.N+ps.M+1):end] = inf_out[(2*ps.N+ps.M+1):end]

  # Loop over metabolites
  for j in 1:(ps.M)
    # Skip metabolites with zero concentration
    if inf_out[ps.N+j] == 0.0
      continue
    end

    # Make vector of +ve and -ve contributions to flux from strains
    pvf = zeros(ps.N)
    nvf = zeros(ps.N)

    # Loop over strains to find values
    for k in 1:(ps.N)
      # Loop over strains reactions
      for l in 1:(ps.mics[k].R)
        # Find reaction
        r = ps.reacs[ps.mics[k].Reacs[l]]
        # Find amount of
        E = Eα(inf_out[2*ps.N+ps.M+k], ps.mics[k], l)
        S_val = [inf_out[ps.N + idx] for idx in r.Rct]
        P_val = [inf_out[ps.N + idx] for idx in r.Prd]
        # Check if reaction has metabolite as a reactant
        if j in r.Rct
          rct_idx = findfirst(==(j), r.Rct)
          q = qs(S_val, P_val, E, l, ps.mics[k], ps.T, r)
          nvf[k] += q * r.stoichR[rct_idx]
        end
        # Or as a product
        if j in r.Prd
          prd_idx = findfirst(==(j), r.Prd)
          q = qs(S_val, P_val, E, l, ps.mics[k], ps.T, r)
          pvf[k] += q * r.stoichP[prd_idx]
        end
      end
    end

    # Calculate approximate perturbation size
    dC = (pvf .- nvf) / ps.δ[j]

    # Loop over metabolites to make perturbation
    for k in 1:(ps.M)
      if k != j
        pout[ps.N+k] = inf_out[ps.N+k]
      else
        pout[ps.N+k] = inf_out[ps.N+k] * (1.0 + frc)
      end
    end

    # Find perturbed forces
    nF = nForce(F, pout, ps)
    # Only interested in forces on ATP
    Fatp[:, j] = nF[(ps.N+ps.M+1):(2*ps.N+ps.M)] .-
                 Fi[(ps.N+ps.M+1):(2*ps.N+ps.M)]

    # Preallocate positive/negative/zero vector
    pn0 = zeros(Int64, ps.N)
    # Preallocate supply/consumption/zero vector
    sc0 = zeros(Int64, ps.N)

    # Loop over strains
    for k in 1:(ps.N)
      # Check whether effect on strain is +ve, -ve or ~0
      # And whether strain has reaction using metabolite j
      if Fatp[k, j] > 1e-5 && j ∈ mts[k]
        # Positive denoted by 1
        pn0[k] = 1
      elseif Fatp[k, j] < -1e-5 && j ∈ mts[k]
        # Negative denoted by -1
        pn0[k] = -1
      end

      # Check net effect of strain
      if dC[k] > 0.0
        # Positive denoted by 1
        sc0[k] = 1
      elseif dC[k] < 0.0
        # Negative denoted by -1
        sc0[k] = -1
      end
    end

    # Loop over vectors to make interaction matrix
    for k in 1:(ps.N)
      for l in 1:(ps.N)
        # Check for positive cases
        if pn0[k] == 1
          # Crossfeeding
          if sc0[l] == 1
            ints[k, l, j] = 2
            in_str[k, l, j] = abs(Fatp[k, j]) * abs(dC[l])
            # Competition
          elseif sc0[l] == -1
            ints[k, l, j] = 1
            in_str[k, l, j] = abs(Fatp[k, j]) * abs(dC[l])
          end
          # And negative case
        elseif pn0[k] == -1
          # Competition via pollution
          if sc0[l] == 1
            ints[k, l, j] = 4
            in_str[k, l, j] = abs(Fatp[k, j]) * abs(dC[l])
            # Obligate syntrophy
          elseif sc0[l] == -1
            ints[k, l, j] = 3
            in_str[k, l, j] = abs(Fatp[k, j]) * abs(dC[l])
          end
        end
      end
      # Rescale all strengths by population, and perturbation size
      in_str[k, :, j] = (inf_out[k] / (ps.mics[k].ρ * ps.mics[k].MC)) *
                        (1 / (frc * inf_out[ps.N+j])) * in_str[k, :, j]
    end
  end

  # Rescale to moles and fraction
  in_str = in_str / (NA)

  # Want to find net reactions, loop over strains
  for j in 1:(ps.N)
    for k in 1:(ps.N)
      # Loop over every metabolite
      for l in 1:(ps.M)
        # Check what type of interaction (if any)
        # effect of k on j
        if ints[j, k, l] == 1
          net_in[j, k] -= in_str[j, k, l]
        elseif ints[j, k, l] == 2
          net_in[j, k] += in_str[j, k, l]
        elseif ints[j, k, l] == 3
          net_in[j, k] += in_str[j, k, l]
        elseif ints[j, k, l] == 4
          net_in[j, k] -= in_str[j, k, l]
        end
      end
    end
  end

  return InteractionResult(ints, in_str, net_in, Fatp, frc)
end

"""
    analyse_ints(result::InteractionResult, ps::FullParameters)

Analyze interaction data and print summary statistics.

# Arguments
- `result::InteractionResult`: Interaction analysis results from quantify_ints()
- `ps::FullParameters`: System parameters

# Returns
- `NamedTuple`: Summary statistics
"""
function analyse_ints(result::InteractionResult, ps::FullParameters)
  ints = result.ints
  in_str = result.in_str
  net_in = result.net_in

  # Count interaction types
  n_competition = count(x -> x == 1, ints)
  n_facilitation = count(x -> x == 2, ints)
  n_syntrophy = count(x -> x == 3, ints)
  n_pollution = count(x -> x == 4, ints)
  n_total = n_competition + n_facilitation + n_syntrophy + n_pollution

  # Calculate mean strengths
  mean_competition = n_competition > 0 ? mean(in_str[findall(x -> x == 1, ints)]) : NaN
  mean_facilitation = n_facilitation > 0 ? mean(in_str[findall(x -> x == 2, ints)]) : NaN
  mean_syntrophy = n_syntrophy > 0 ? mean(in_str[findall(x -> x == 3, ints)]) : NaN
  mean_pollution = n_pollution > 0 ? mean(in_str[findall(x -> x == 4, ints)]) : NaN

  # Print summary
  println("\n--- Interaction Analysis Summary ---")
  println("Total interactions: $n_total")
  println("  Competition: $n_competition ($(round(n_competition/n_total*100, digits=1))%)")
  println("  Facilitation: $n_facilitation ($(round(n_facilitation/n_total*100, digits=1))%)")
  println("  Syntrophy: $n_syntrophy ($(round(n_syntrophy/n_total*100, digits=1))%)")
  println("  Pollution: $n_pollution ($(round(n_pollution/n_total*100, digits=1))%)")
  println("\nMean interaction strengths:")
  if !isnan(mean_competition)
    println("  Competition: $(round(mean_competition, sigdigits=3))")
  end
  if !isnan(mean_facilitation)
    println("  Facilitation: $(round(mean_facilitation, sigdigits=3))")
  end
  if !isnan(mean_syntrophy)
    println("  Syntrophy: $(round(mean_syntrophy, sigdigits=3))")
  end
  if !isnan(mean_pollution)
    println("  Pollution: $(round(mean_pollution, sigdigits=3))")
  end
  println("--------------------------------\n")

  return (
    total=n_total,
    competition=n_competition,
    facilitation=n_facilitation,
    syntrophy=n_syntrophy,
    pollution=n_pollution,
    mean_competition=mean_competition,
    mean_facilitation=mean_facilitation,
    mean_syntrophy=mean_syntrophy,
    mean_pollution=mean_pollution
  )
end

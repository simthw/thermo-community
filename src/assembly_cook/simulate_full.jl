# A script to run the dynamics for the full model.
export full_simulate, θ, qs, Eα

export full_simulate_syn

# Buffered metabolites that are treated as infinite reservoirs (constant concentration)
# Index 5 = H2O (water) - effectively unlimited in aqueous systems
# Index 8 = H+ (hydrogen ion) - buffered by carbonate system to maintain stable pH
const BUFFERED_METABOLITES = Set([5, 8])

# function to find the thermodynamic term θ, for reactions with multiple substrates and products
# Returns the thermodynamic inhibition factor [0, 1], where:
#   - θ = 0: no thermodynamic inhibition (reaction proceeds at maximal rate)
#   - θ = 1: complete inhibition (reaction cannot proceed forward)

# θ is computed as Q/Keq and clamped to [0, 1]:
#   - Q = reaction quotient = ∏(P^ν) / ∏(S^ν)
#   - Keq = equilibrium constant = exp((-ΔG0 - η*ΔGATP) / RT)

# Special cases:
#   - When any substrate S <= 0: θ = 1 (no reaction possible without substrate)
#   - When any product P <= 0: θ = 0 (no thermodynamic inhibition, maximal rate)
#   - When Q >= Keq: θ = 1 (reaction at equilibrium or thermodynamically reversed)
#   - When Q << Keq: θ ≈ 0 (far from equilibrium, minimal inhibition)

# Note: θs = Q/Keq can exceed 1 when the reaction is at or past equilibrium.
#       We clamp to 1.0 because the rate formula uses (1 - θ), and negative
#       rates are unphysical in this model (would imply spontaneous reversal).
#       Negative rates are prevented by max(q, 0) in the rate calculation.
function θ(S::Vector{Float64}, P::Vector{Float64}, T::Float64, η::Float64, r::Reaction)
  # Catch perverse cases that sometimes arise
  if any(S .<= 0.0)
    θs = 1.0  # No reaction possible without substrate
  elseif any(P .<= 0.0)
    θs = 0.0  # No thermodynamic inhibition when no product (matches Assembly)
  else
    Q = prod(P .^ r.stoichP) / prod(S .^ r.stoichR)
    θs = Q / Keq(T, η, r.ΔG0)
  end
  # Clamp θ to [0, 1]: values > 1 indicate reaction is at/past equilibrium
  return clamp(θs, 0.0, 1.0)
end

# Extra function to find rate from thermodynamic inhibition, just useful for plotting
function qs(mic::MicrobeP, S::Vector{Float64}, P::Vector{Float64}, E::Float64, θs::Float64)
  limiting = minimum(S)
  q = mic.kc[1] * E * limiting * (1 - θs) / (mic.KS[1] + limiting * (1 + mic.kr[1] * θs))
  return (max(q, 0.0))
end

# function to find the rate of substrate consumption by a particular reaction
function qs(S::Vector{Float64}, P::Vector{Float64}, E::Float64, i::Int64, ps::MicrobeP, T::Float64,
  r::Reaction; thermodynamic_feedback::Bool=true)
  θs = thermodynamic_feedback ? θ(S, P, T, ps.η[i], r) : 0.0
  limiting = minimum(S ./ r.stoichR)
  q = ps.kc[i] * E * limiting * (1 - θs) / (ps.KS[i] + limiting * (1 + ps.kr[i] * θs))
  # Ensure that negative value cannot be returned
  return (max(q, 0.0))
end

# function to calculate the amount of a particular enzyme a strain has
function Eα(ϕR::Float64, ps::MicrobeP, i::Int64)
  E = ps.MC * max(1 - ϕR - ps.ϕH, 0.0) * ps.ϕP[i] / (ps.n[2])
  return (E)
end

# function to find the elongation rate γ
function γs(a::Float64, ps::MicrobeP)
  γ = ps.γm * a / (ps.Kγ + a)
  return (γ)
end

# function to find the growth rate λ
function λs(a::Float64, ϕR::Float64, ps::MicrobeP)
  γ = γs(a, ps)
  λ = (γ * ϕR * ps.Pb) / ps.n[1]
  return (λ)
end

# function to find ϕR based on the energy concentration
function ϕ_R(a::Float64, ps::MicrobeP)
  ϕ = (1 - ps.ϕH) * a / (ps.KΩ + a)
  return (ϕ)
end

# function to implement the consumer resource dynamics
function full_dynamics!(dx::Array{Float64,1}, x::Array{Float64,1}, ps::FullParameters,
  rate::Array{Float64,2}, t::Float64, thermodynamic_feedback::Bool=true)
  fill!(dx, 0.0)
  fill!(rate, 0.0)

  for i in 1:ps.N
    energy = x[ps.N+ps.M+i]
    if x[i] <= 1e-10 || energy <= 0.33
      continue
    end
    mic = ps.mics[i]
    for k in 1:mic.R
      reaction_id = mic.Reacs[k]
      reaction = ps.reacs[reaction_id]
      enzyme = Eα(x[2*ps.N+ps.M+i], mic, k)
      substrates = [max(x[ps.N+rct], 0.0) for rct in reaction.Rct]
      products = [max(x[ps.N+prd], 0.0) for prd in reaction.Prd]
      rate[i, reaction_id] = qs(
        substrates,
        products,
        enzyme,
        k,
        mic,
        ps.T,
        reaction;
        thermodynamic_feedback=thermodynamic_feedback,
      )
    end
  end

  for i in 1:(ps.N)
    if x[i] <= 1e-10 || x[ps.N+ps.M+i] <= 0.33
      dx[i] = 0.0
      dx[ps.N+ps.M+i] = 0.0
      dx[2*ps.N+ps.M+i] = 0.0
    else
      phi = max(x[2*ps.N+ps.M+i], 0.0)
      λ = λs(x[ps.N+ps.M+i], phi, ps.mics[i])
      dx[i] = (λ - ps.mics[i].d) * x[i]
      ϕR = ϕ_R(x[ps.N+ps.M+i], ps.mics[i])
      τ = ps.mics[i].fd / max(λ, eps(Float64))
      dx[2*ps.N+ps.M+i] = (ϕR - x[2*ps.N+ps.M+i]) / τ
      # Calculate energy production rate J = Σ(η_j · rate_j)
      # All reactions (including endergonic H2 production) use η·rate formulation
      #
      # Endergonic reaction handling (ΔG0 > 0, e.g., syntrophic H2 production):
      # - Default η = 0 for standalone endergonic reactions
      # - For H2 production reactions coupled with H2 consumers (methanogens/sulfate reducers),
      #   a conservative positive η is assigned if the combined system is exergonic
      # - Thermodynamic feasibility is always controlled by θ via Keq
      # - This allows syntrophic organisms to gain energy from H2 production when
      #   partnered with H2 consumers, while preventing unrealistic energy generation
      J = 0.0
      for j in 1:(ps.mics[i].R)
        rj = ps.mics[i].Reacs[j]
        q_j = rate[i, rj]
        J += ps.mics[i].η[j] * q_j
      end
      dx[ps.N+ps.M+i] = J - (ps.mics[i].MC * ps.mics[i].ρ + x[ps.N+ps.M+i]) * λ
    end
  end
  for i in (ps.N+1):(ps.N+ps.M)
    met_idx = i - ps.N
    # Buffered metabolites (e.g., water) have constant concentration - no dilution
    if met_idx ∈ BUFFERED_METABOLITES
      dx[i] = 0.0
    else
      dx[i] = ps.κ[met_idx] - ps.δ[met_idx] * max(x[i], 0.0)
    end
  end

  for i in 1:(ps.N)
    for j in 1:(ps.mics[i].R)
      r = ps.reacs[ps.mics[i].Reacs[j]]
      for (k, prd) in enumerate(r.Prd)
        # Skip buffered metabolites - their concentration is maintained externally
        if prd ∉ BUFFERED_METABOLITES
          dx[ps.N+prd] += rate[i, ps.mics[i].Reacs[j]] * r.stoichP[k] * max(x[i], 0.0) / NA
        end
      end
      for (k, rct) in enumerate(r.Rct)
        # Skip buffered metabolites - their concentration is maintained externally
        if rct ∉ BUFFERED_METABOLITES
          dx[ps.N+rct] -= rate[i, ps.mics[i].Reacs[j]] * r.stoichR[k] * max(x[i], 0.0) / NA
        end
      end
    end
  end
  # Keep the vector field inside the non-negative domain without mutating the
  # solver-owned state. Positive production remains possible from a zero pool.
  for i in (ps.N+1):(ps.N+ps.M)
    if x[i] <= 0.0 && dx[i] < 0.0
      dx[i] = 0.0
    end
  end
  for i in 1:ps.N
    if x[i] <= 0.0 && dx[i] < 0.0
      dx[i] = 0.0
    end
  end
  for i in (2*ps.N+ps.M+1):(3*ps.N+ps.M)
    if x[i] <= 0.0 && dx[i] < 0.0
      dx[i] = 0.0
    end
  end
  return (dx)
end

# Simulation code to run one instance of the simulation with a user defined starting condition
# verbose: if true, print solver diagnostics (default: true)
function full_simulate(ps::FullParameters, Tmax::Float64, pop::Array{Float64,1},
  conc::Array{Float64,1}, as::Array{Float64,1}, ϕs::Array{Float64,1};
  verbose::Bool=true, thermodynamic_feedback::Bool=true, return_metadata::Bool=false)
  @assert length(pop) == ps.N "From parameter set expected $(ps.N) strains"
  @assert length(conc) == ps.M "From parameter set expected $(ps.M) metabolites"
  @assert length(as) == length(pop) "Every strain must have an energy concentration"
  @assert length(ϕs) == length(pop) "Every strain must have a ribosome fraction"
  rate = zeros(ps.N, ps.O)
  dyns!(dx, x, ps, t) = full_dynamics!(dx, x, ps, rate, t, thermodynamic_feedback)
  tspan = (0, Tmax)
  x0 = [pop; conc; as; ϕs]
  prob = ODEProblem(dyns!, x0, tspan, ps)
  sol = DifferentialEquations.solve(
    prob, Rodas5(autodiff=AutoFiniteDiff()); reltol=1e-6, abstol=1e-8, maxiters=1e7)
  if verbose
    println(sol.destats)
    println(sol.retcode)
  end
  metadata = (
    retcode=string(sol.retcode),
    successful=SciMLBase.successful_retcode(sol.retcode),
    final_time=Float64(sol.t[end]),
    reached_tmax=isapprox(Float64(sol.t[end]), Tmax; rtol=1e-8, atol=1e-8),
    thermodynamic_feedback=thermodynamic_feedback,
  )
  result = (Array(sol'), sol.t)
  return return_metadata ? (result..., metadata) : result
end

# Same as the above but with a specific strain given a massively increased death rate
# This should help us investigate syntrophic pairs by adding array 'inds'
function full_simulate_syn(ps::FullParameters, Tmax::Float64, pop::Array{Float64,1},
  conc::Array{Float64,1}, as::Array{Float64,1}, ϕs::Array{Float64,1}, inds::Array{Int64,1};
  thermodynamic_feedback::Bool=true)
  @assert length(pop) == ps.N "From parameter set expected $(ps.N) strains"
  @assert length(conc) == ps.M "From parameter set expected $(ps.M) metabolites"
  @assert length(as) == length(pop) "Every strain must have an energy concentration"
  @assert length(ϕs) == length(pop) "Every strain must have a ribosome fraction"
  rate = zeros(ps.N, ps.O)
  nmics = Array{MicrobeP,1}(undef, length(ps.mics))
  for i in 1:length(ps.mics)
    if i ∉ inds
      nmics[i] = ps.mics[i]
    else
      m = ps.mics[i]
      nd = 25.0 * m.d
      nmics[i] = make_MicrobeP(
        m.MC, m.γm, m.ρ, m.Kγ, m.Pb, nd, m.ϕH, m.KΩ, m.fd, m.R,
        m.Reacs, m.η, m.kc, m.KS, m.kr, m.n, m.ϕP)
    end
  end
  ps2 = make_FullParameters(ps.N, ps.M, ps.O, ps.T, ps.κ, ps.δ, ps.reacs, nmics)
  dyns!(dx, x, ps2, t) = full_dynamics!(dx, x, ps2, rate, t, thermodynamic_feedback)
  tspan = (0, Tmax)
  x0 = [pop; conc; as; ϕs]
  prob = ODEProblem(dyns!, x0, tspan, ps2)
  sol = DifferentialEquations.solve(
    prob, Rodas5(autodiff=AutoFiniteDiff()); reltol=1e-6, abstol=1e-8, maxiters=1e7)
  return (Array(sol'), sol.t)
end

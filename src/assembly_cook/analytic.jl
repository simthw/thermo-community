# Script that provides analytic functions such as Lyapunov exponents for use in simulations

using SymPy: @syms, subs, diff, Sym

export Force, nForce, Jacobian, nJacobian

# function to return value of theta between 0 and 1
function bound_θ(S::Vector{Float64}, P::Vector{Float64}, T::Float64, η::Float64, r::Reaction)
  if any(S .<= 0.0)
    θs = 1.0
  elseif any(P .<= 0.0)
    θs = 0.0
  else
    θs = prod(P .^ r.stoichP) / prod(S .^ r.stoichR)
    θs /= Keq(T, η, r.ΔG0)
  end
  return min(1.0, θs)
end

# function to return rate in symbolic form
function symb_rate(S, θ, stoichR::Vector{Float64})
  # Define kinetic/thermodynamic parameters
  kc, KS, kr = symbols("kc, KS, kr")
  limiting = if length(S) == 1
    S[1] / stoichR[1]
  else
    SymPy.sympy.Min((S[i] / stoichR[i] for i in eachindex(S))...)
  end
  # Start with numerator
  n = kc * limiting * (1 - θ)
  # Then find denominator
  d = KS + limiting * (1 + kr * θ)
  # Divide numerator by denominator
  q = n / d
  return (q)
end

theta_name(rc::Reaction, strain_idx::Int64) = "θR$(rc.ID)N$(strain_idx)"
theta_symbol(rc::Reaction, strain_idx::Int64) = symbols(theta_name(rc, strain_idx))

function reaction_substrate_symbols(rc::Reaction)
  return [symbols("M$(idx)") for idx in rc.Rct]
end

function reaction_product_symbols(rc::Reaction)
  return [symbols("M$(idx)") for idx in rc.Prd]
end

function symbolic_bound_θ(rc::Reaction, η::Float64, T::Float64)
  S = reaction_substrate_symbols(rc)
  P = reaction_product_symbols(rc)
  K = Keq(T, η, rc.ΔG0)
  θ_expr = prod(P[idx]^rc.stoichP[idx] for idx in eachindex(P)) /
           (prod(S[idx]^rc.stoichR[idx] for idx in eachindex(S)) * K)
  return SymPy.sympy.Min(θ_expr, 1.0)
end

# function to find forces symbolically
function Force(ps::FullParameters, F::Array{Sym,1})
  # Check Force provided is right size
  @assert length(F) == 3 * ps.N + ps.M "Preallocated force vector incorrect size"
  # Define ecosystem and microbe level symbols, create symbolic variables
  N, η, κ, M, δ, d, γ, a, Kγ, Pb, ϕR, nr, np = symbols("N, η, κ, M, δ, d, γ, a, Kγ, fb, ϕR, nr, np")
  MC, ϕH, ρ, KΩ, fd = symbols("MC, ϕH, ρ, KΩ, fd")
  # And reaction level symbols
  kc, KS, kr, ϕP = symbols("kc, KS, kr, ϕP")
  # Loop over all microbes to begin with, dNi/dt
  for i in 1:(ps.N)
    ex = -d # death rate (symbolic)
    ex += (γ * a / (a + Kγ)) * ϕR * Pb / nr # add growth rate (symbolic)
    ex *= N # multiply by population size (symbolic)
    # Sub in strain/environment level variables
    ex = subs(ex, γ => ps.mics[i].γm, Kγ => ps.mics[i].Kγ, Pb => ps.mics[i].Pb,
      nr => ps.mics[i].n[1])
    F[i] = subs(ex, d => ps.mics[i].d, a => "a$(i)", N => "N$(i)", ϕR => "ϕR$(i)")
  end
  # Then loop over metabolite, dMi/dt = influx - decay - consumption + production
  for i in (ps.N+1):(ps.N+ps.M)
    # Set metabolite number, Mn
    Mn = i - ps.N # from 1 to ps.M
    # Add external dynamics
    ex = κ - δ * M
    # Sub in relevant values
    ex = subs(ex, κ => ps.κ[Mn], δ => ps.δ[Mn], M => "M$(Mn)")
    # Loop over all strains, strain loop j, loops over each microbial strain (1 to ps.N), and reaction loop k, 
    # for each strain j, loops over its reactions (1 to ps.mics[j].R, where R is number of reactions strain j can perform)
    for j in 1:(ps.N)
      # Then loop over each reaction for each strain
      for k in 1:(ps.mics[j].R)
        rc = ps.reacs[ps.mics[j].Reacs[k]]
        # Check if metabolite is used as reactant, Case 1: Metabolite is Reactant (nR == Mn)
        rct_idx = findfirst(==(Mn), rc.Rct)
        if !isnothing(rct_idx)
          S = reaction_substrate_symbols(rc)
          θ = theta_symbol(rc, j)
          # Then make symbolic rate
          q = symb_rate(S, θ, rc.stoichR)
          # Then multiply by population, and amount of enzyme
          ex -= rc.stoichR[rct_idx] * q * N * MC * ϕP * (1 - ϕR - ϕH) / (np * NA)
          # sub in microbe identifier and parameters
          ex = subs(ex, N => "N$(j)", MC => ps.mics[j].MC, ϕR => "ϕR$(j)")
          ex = subs(ex, ϕH => ps.mics[j].ϕH, np => ps.mics[j].n[2])
          # Then sub in all metabolite level parameters
          ex = subs(ex, kc => ps.mics[j].kc[k], KS => ps.mics[j].KS[k])
          ex = subs(ex, kr => ps.mics[j].kr[k], ϕP => ps.mics[j].ϕP[k])
          # Or if metabolite is produced as product
        end

        prd_idx = findfirst(==(Mn), rc.Prd)
        if !isnothing(prd_idx)
          S = reaction_substrate_symbols(rc)
          θ = theta_symbol(rc, j)
          # Then make symbolic rate
          q = symb_rate(S, θ, rc.stoichR)
          # Then multiply by population, and amount of enzyme
          ex += rc.stoichP[prd_idx] * q * N * MC * ϕP * (1 - ϕR - ϕH) / (np * NA)
          # sub in microbe identifier and parameters
          ex = subs(ex, N => "N$(j)", MC => ps.mics[j].MC, ϕR => "ϕR$(j)")
          ex = subs(ex, ϕH => ps.mics[j].ϕH, np => ps.mics[j].n[2])
          # Then sub in all metabolite level parameters
          ex = subs(ex, kc => ps.mics[j].kc[k], KS => ps.mics[j].KS[k])
          ex = subs(ex, kr => ps.mics[j].kr[k], ϕP => ps.mics[j].ϕP[k])
        end
      end
    end
    # No further substitutions required
    F[i] = ex
  end
  # Next loop over energy concentrations
  for i in (ps.N+ps.M+1):(2*ps.N+ps.M)
    # Find relevant microbe
    mic = ps.mics[i-ps.N-ps.M] # from 1 to ps.N
    # energy consumption rate due to growth, negative sign
    ex = -(γ * a / (a + Kγ)) * Pb * ϕR / nr
    # Multiply it by translation cost and concentration
    ex *= (ρ * MC + a)
    # Then loop over reactions
    for j in 1:(mic.R)
      # Find jth reaction
      rc = ps.reacs[mic.Reacs[j]]
      S = reaction_substrate_symbols(rc)
      θ = theta_symbol(rc, i - ps.N - ps.M)
      # Then make symbolic rate
      q = symb_rate(S, θ, rc.stoichR)
      # Multiple expression by relevant factors, including amount of enzyme
      ex += η * q * MC * ϕP * (1 - ϕR - ϕH) / (np)
      # Now sub in all of reaction level parameters
      ex = subs(ex, kc => mic.kc[j], KS => mic.KS[j], kr => mic.kr[j], η => mic.η[j],
        ϕP => mic.ϕP[j])
    end
    # Sub in all of relevant variables
    ex = subs(ex, γ => mic.γm, a => "a$(i-ps.N-ps.M)", Kγ => mic.Kγ, Pb => mic.Pb,
      nr => mic.n[1])
    F[i] = subs(ex, ϕR => "ϕR$(i-ps.N-ps.M)", ρ => mic.ρ, MC => mic.MC, ϕH => mic.ϕH,
      np => mic.n[2])
  end
  # Finally do same for ribosome fractions
  for i in (2*ps.N+ps.M+1):(3*ps.N+ps.M)
    # Find relevant microbe
    mic = ps.mics[i-2*ps.N-ps.M]
    # Subtract current ribosome fraction
    ex = -ϕR
    # Add "optimal" ribosome fraction
    ex += (a / (a + KΩ)) * (1 - ϕH)
    # Then divide by characteristic time scale
    ex /= fd * nr / ((γ * a / (a + Kγ)) * Pb * ϕR)
    # Sub in all of relevant variable
    ex = subs(ex, ϕR => "ϕR$(i-2*ps.N-ps.M)", a => "a$(i-2*ps.N-ps.M)", KΩ => mic.KΩ)
    F[i] = subs(ex, ϕH => mic.ϕH, fd => mic.fd, nr => mic.n[1], γ => mic.γm,
      Kγ => mic.Kγ, Pb => mic.Pb)
  end
  return (F)
end

# function to find numerical values of forces at a particular point in  space
function nForce(F::Array{Sym,1}, C::Array{Float64,1}, ps::FullParameters)
  # Copy F to a new object to prevent overwriting
  f = copy(F)
  # Sub in steady state population values
  for i in 1:(ps.N)
    for j in 1:(3*ps.N+ps.M)
      f[j] = subs(f[j], "N$i" => C[i])
    end
  end
  # Sub in steady state metabilites concentrations values
  for i in (ps.N+1):(ps.N+ps.M)
    for j in 1:(3*ps.N+ps.M)
      f[j] = subs(f[j], "M$(i-ps.N)" => C[i])
    end
  end
  # Then sub in theta values
  for i in 1:(ps.N)
    for j in 1:(ps.mics[i].R)
      rc = ps.reacs[ps.mics[i].Reacs[j]]
      S = [C[ps.N + idx] for idx in rc.Rct]
      P = [C[ps.N + idx] for idx in rc.Prd]
      θ = bound_θ(S, P, ps.T, ps.mics[i].η[j], rc)
      for k in 1:(3*ps.N+ps.M)
        f[k] = subs(f[k], theta_name(rc, i) => θ)
      end
    end
  end
  # Next sub in energies
  for i in (ps.N+ps.M+1):(2*ps.N+ps.M)
    for j in 1:(3*ps.N+ps.M)
      f[j] = subs(f[j], "a$(i-ps.N-ps.M)" => C[i])
    end
  end
  # Finally sub in ribosome fractions
  for i in (2*ps.N+ps.M+1):(3*ps.N+ps.M)
    for j in 1:(3*ps.N+ps.M)
      f[j] = subs(f[j], "ϕR$(i-2*ps.N-ps.M)" => C[i])
    end
  end
  # convert vector into a float
  f = convert(Array{Float64}, f)
  return (f)
end

# function to find an analytic form of the Jacobian
function Jacobian(F::Array{Sym,1}, J::Array{Sym,2}, ps::FullParameters,
  cons::Array{Float64,1})
  # Check Force provided is right size
  @assert size(J) == (length(F), length(F)) "preallocated Jacobian wrong size"
  @assert size(J) == (3 * ps.N + ps.M, 3 * ps.N + ps.M) "preallocated Jacobian wrong size"
  @assert length(cons) == ps.M "concentrations wrong length"

  # Copy F to a new object to prevent overwriting
  f = copy(F)
  # Loop over strains to sub in theta values
  for i in 1:(ps.N)
    # For each microbe loop over reactions it has
    for j in 1:(ps.mics[i].R)
      # Find reaction
      rc = ps.reacs[ps.mics[i].Reacs[j]]
      # And corresponding θ symbol
      θ = theta_name(rc, i)
      # Check if substrate concentration is zero
      if any(cons[idx] <= 0.0 for idx in rc.Rct)
        # In this case substitute 1 in to remove  relevant terms
        for k in 1:(3*ps.N+ps.M)
          f[k] = subs(f[k], θ => 1.0)
        end
      elseif any(cons[idx] <= 0.0 for idx in rc.Prd)
        for k in 1:(3*ps.N+ps.M)
          f[k] = subs(f[k], θ => 0.0)
        end
      else
        θ_expr = symbolic_bound_θ(rc, ps.mics[i].η[j], ps.T)
        # Sub true theta value into  vector of forces
        for k in 1:(3*ps.N+ps.M)
          f[k] = subs(f[k], θ => θ_expr)
        end
      end
    end
  end
  # Now loop over variables to differentiate by
  for j in 1:(3*ps.N+ps.M)
    # Find variable to differentiate by
    if j <= ps.N
      dx = symbols("N$j")
    elseif j <= ps.N + ps.M
      dx = symbols("M$(j-ps.N)")
    elseif j <= 2 * ps.N + ps.M
      dx = symbols("a$(j-ps.N-ps.M)")
    else
      dx = symbols("ϕR$(j-2*ps.N-ps.M)")
    end
    for i in 1:(3*ps.N+ps.M)
      J[i, j] = diff(f[i], dx)
    end
  end
  return (J)
end

# function to find numerical values of forces at a particular point in  space
function nJacobian(J::Array{Sym,2}, C::Array{Float64,1}, ps::FullParameters)
  # Copy F to a new object to prevent overwriting
  Jc = copy(J)
  # Sub in steady state population values
  for i in 1:(ps.N)
    for j in 1:(3*ps.N+ps.M)
      for k in 1:(3*ps.N+ps.M)
        Jc[j, k] = subs(Jc[j, k], "N$i" => C[i])
      end
    end
  end
  # Sub in steady state concentrations values
  for i in (ps.N+1):(ps.N+ps.M)
    for j in 1:(3*ps.N+ps.M)
      for k in 1:(3*ps.N+ps.M)
        Jc[j, k] = subs(Jc[j, k], "M$(i-ps.N)" => C[i])
      end
    end
  end
  # Next sub in energies
  for i in (ps.N+ps.M+1):(2*ps.N+ps.M)
    for j in 1:(3*ps.N+ps.M)
      for k in 1:(3*ps.N+ps.M)
        Jc[j, k] = subs(Jc[j, k], "a$(i-ps.N-ps.M)" => C[i])
      end
    end
  end
  # Finally sub in ribosome fractions
  for i in (2*ps.N+ps.M+1):(3*ps.N+ps.M)
    for j in 1:(3*ps.N+ps.M)
      for k in 1:(3*ps.N+ps.M)
        Jc[j, k] = subs(Jc[j, k], "ϕR$(i-2*ps.N-ps.M)" => C[i])
      end
    end
  end
  # convert vector into a float
  Jc = convert(Array{Float64}, Jc)
  return (Jc)
end

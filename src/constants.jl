export Rgas, ΔGATP, NA
export mvector, Keq

using Distributions: Normal

const Rgas = 8.31446261815324
const ΔGATP = 75_000.0
const NA = 6.02214076e23

function mvector(N::Int64, mm::Float64, sdm::Float64)
  @assert mm - 5*sdm >= 0.0 "Negative maintenance energies are not allowed"
  return rand(Normal(mm, sdm), N)
end

Keq(T::Float64, η::Float64, ΔG0::Float64) = exp((-ΔG0 - η * ΔGATP) / (Rgas * T))

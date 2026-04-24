using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

function evolve_deepcliffords(ψ::CAMPS,
                              t::Integer,
                              paulistrings::Vector{<:PauliOperator},
                              phases::Vector{<:Real};
                              showprogress = false)
  k = 0
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    applyQOp!(ψ, QOp(:randCliffCircuit))

    P = paulistrings[s]
    ϕ = phases[s]
    N = length(ψ)
    I = PauliOperator(0x0, fill(false,N), fill(false,N))
    C = inv(ψ.Cdag)

    nature = paulinature(k, C, P)
    if nature == :disentanglable
      ψ.Cdag *= inv(disentangler(k, C, P))
      addmagicstate!(ψ, k, ϕ)
      k += 1
    elseif nature == :logical
      R = PauliSum([cos(ϕ), sin(ϕ)], Stabilizer([I,P]))
      applyGate!(ψ, R)
    elseif nature == :trivial
      next!(progressbar)
      continue
    end
    next!(progressbar)
  end
  return ψ, k
end

"Turn the (k+1)th qubit from |0⟩ to the Liu and Clark (2025) magic state"
function addmagicstate!(ψ::CAMPS, k::Integer, phase::Real)
  magifier_os = OpSum()
  magifier_os += cos(phase), "Id", k+1
  magifier_os += -im * sin(phase), "X", k+1
  sites = siteinds(ψ.mps)
  magifier = MPO(magifier_os, sites)
  ψ.mps = ITensors.apply(magifier, ψ.mps)
  return nothing
end
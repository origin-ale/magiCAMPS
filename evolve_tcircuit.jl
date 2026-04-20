using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

function evolve_tcircuit(ψ::CAMPS, t::Int; showprogress = false)
  k = 0
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    applyQOp!(ψ, QOp(:randCliffCircuit))

    C = inv(ψ.Cdag)
    nature = paulinature(k, C, Z₁)
    if nature == :disentanglable
      ψ.Cdag *= inv(disentangler(k, C, Z₁))
      addmagicstate!(ψ, k)
      k += 1
    elseif nature == :logical
      applyGate!(ψ, T₁)
    elseif nature == :trivial
      next!(progressbar)
      continue
    end
    next!(progressbar)
  end
  return ψ, k
end

"Turn the (k+1)th qubit from |0⟩ to the Liu and Clark (2025) magic state"
function addmagicstate!(ψ::CAMPS, k::Integer)
  magifier_os = OpSum()
  magifier_os += cos(π/8), "Id", k+1
  magifier_os += -im * sin(π/8), "X", k+1
  sites = siteinds(ψ.mps)
  magifier = MPO(magifier_os, sites)
  ψ.mps = ITensors.apply(magifier, ψ.mps)
  return nothing
end
using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

function evolve_tcircuit(ψ::CAMPS, t::Int)
  k = 0
  @showprogress desc = "Evolving…" for s in 1:t
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
      continue
    end
  end
  return ψ, k
end

function addmagicstate!(ψ::CAMPS, k::Integer)
  magifier_os = OpSum()
  magifier_os += cos(π/8), "Id", k+1
  magifier_os += -im * sin(π/8), "X", k+1
  sites = siteinds(ψ.mps)
  magifier = MPO(magifier_os, sites)
  ψ.mps = ITensors.apply(magifier, ψ.mps)
  return nothing
end
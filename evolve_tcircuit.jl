using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

function evolve_tcircuit(ψ::CAMPS, t::Int)
  k = 0
  @showprogress desc = "Evolving…" for s in 1:t # At each timestep…
    applyQOp!(ψ, QOp(:randCliffCircuit)) # Apply the deep Clifford circuit

    C = inv(ψ.Cdag)
    nature = paulinature(k, C, Z₁)
    if nature == :disentanglable
      ψ.Cdag *= inv(disentangler(k, C, Z₁))
      k += 1
    elseif nature == :logical
      applyGate!(ψ, T₁)
    elseif nature == :trivial
      continue
    end
  end
  return ψ, k
end
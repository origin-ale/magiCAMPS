include("algorithm.jl")

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

N = 10
xbits = fill(false, N)
zbits = fill(false, N)
zbits[1] = true
Z₁ = PauliOperator(0x0, xbits, zbits)
T₁ = PauliSum(QOp(:T), N, 1)
t = 10
ψ = CAMPS(N) # Initialize CAMPS
ψ_evo, k = evolve_tcircuit(ψ, t)
println("$(N-k) free qubit(s) left.")
println(ψ_evo)

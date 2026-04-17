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

N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
t = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10
xbits = fill(false, N)
zbits = fill(false, N)
zbits[1] = true
Z₁ = PauliOperator(0x0, xbits, zbits)
T₁ = PauliSum(QOp(:T), N, 1)
ψ = CAMPS(N) # Initialize CAMPS
ψ_evo, k = evolve_tcircuit(ψ, t)

println("$(N-k) free qubit(s) left.")
print("Final MPS:")
println(ψ_evo.mps)
avg_ee = maximum(eEntropys!(ψ_evo.mps))/N
println("Avg entanglement entropy $(avg_ee)")
avg_sre = sEntropy(ψ_evo.mps, N^2)
println("Avg 2-SRE $(avg_sre)")

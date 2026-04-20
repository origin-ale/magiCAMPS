include("algorithm.jl")
include("evolve_tcircuit.jl")

using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
t = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10
xbits = fill(false, N)
zbits = fill(false, N)
zbits[1] = true
Z₁ = PauliOperator(0x0, xbits, zbits)
T₁ = PauliSum(QOp(:T), N, 1)

ψ = CAMPS(N) # Initialize CAMPS
ψ_evo, k = evolve_tcircuit(ψ, t; showprogress = true)

println("$(N-k) free qubit(s) left.")
print("Final MPS:")
println(ψ_evo.mps)
avg_ee = maximum(eEntropys!(ψ_evo.mps))/N
println("Avg entanglement entropy $(avg_ee)")
avg_sre = sEntropy(ψ_evo.mps, N^2; α =2)/N
println("Avg 2-SRE $(avg_sre)")

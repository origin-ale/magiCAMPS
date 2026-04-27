using DisentangleCAMPS
using CliffordMPS

using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
t = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10
χ = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 64

ψ = CAMPS(N)
paulistrings = [PauliOperator(0x0, rand(Bool,N), rand(Bool,N)) for _ in 1:t]
phases = 2π*rand(Float64, (t,))

ψ_evo, k, tstop = DisentangleCAMPS.evolve_bonddim(ψ, χ, paulistrings, phases; showprogress = true)

println("Evolution stopped at t=$(tstop).")
println("$(N-k) free qubit(s) left.")
print("Final MPS:")
println(ψ_evo.mps)
avg_ee = maximum(eEntropys!(ψ_evo.mps))/N
println("Avg entanglement entropy $(avg_ee)")
avg_sre = sEntropy(ψ_evo.mps, N^2; α =2)/N
println("Avg 2-SRE $(avg_sre)")
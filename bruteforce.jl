using QuantumClifford
using LinearAlgebra
using Random

# Brute-force statevector simulation of Pauli rotations.

function pauli_to_matrix(P::PauliOperator)
  N = length(P)
  σI = ComplexF64[1 0; 0 1]
  σX = ComplexF64[0 1; 1 0]
  σY = ComplexF64[0 -im; im 0]
  σZ = ComplexF64[1 0; 0 -1]
  M = ComplexF64[1.0;;]
  for j in 1:N
    x, z = P[j]
    σ = (x && z) ? σY : (x ? σX : (z ? σZ : σI))
    M = kron(M, σ)
  end
  phase_factor = (-1)^(P.phase[1] >> 1) * (im)^(P.phase[1] & 1)
  return phase_factor * M
end

function rotate_state(ψ::Vector{ComplexF64}, P::PauliOperator, ϕ::Real)
  M = pauli_to_matrix(P)
  return cos(ϕ/2) * ψ - im * sin(ϕ/2) * (M * ψ)
end

function expect_pauli(ψ::Vector{ComplexF64}, P::PauliOperator)
  M = pauli_to_matrix(P)
  return real(ψ' * M * ψ)
end

Random.seed!(42)
N = 5
T = 5
paulistrings = [PauliOperator(0x0, rand(Bool,N), rand(Bool,N)) for _ in 1:T]
phases = 2π*rand(Float64, T)

function go()
  Random.seed!(42)
  N = 5
  paulistrings = [PauliOperator(0x0, rand(Bool,N), rand(Bool,N)) for _ in 1:5]
  phases = 2π*rand(Float64, 5)
  ψ = zeros(ComplexF64, 2^N)
  ψ[1] = 1.0
  for s in 1:3
    ψ = rotate_state(ψ, paulistrings[s], phases[s])
  println("=== Step $s P=$(paulistrings[s]) ϕ=$(phases[s]) ===")
  for j in 1:N
    pX = PauliOperator(0x0, [i==j for i in 1:N], fill(false,N))
    pY = PauliOperator(0x0, [i==j for i in 1:N], [i==j for i in 1:N])
    pZ = PauliOperator(0x0, fill(false,N), [i==j for i in 1:N])
    ex = round(expect_pauli(ψ, pX), digits=4)
    ey = round(expect_pauli(ψ, pY), digits=4)
    ez = round(expect_pauli(ψ, pZ), digits=4)
    println("  qubit $j: ⟨X⟩=$ex ⟨Y⟩=$ey ⟨Z⟩=$ez")
  end
  end
end

go()

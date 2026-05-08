using DisentangleCAMPS
using DisentangleCAMPS: paulinature, disentangler, apply
using CliffordMPS
using QuantumClifford
using Random

function go(mode)
  Random.seed!(42)
  N = 5
  T = 5
  paulistrings = [PauliOperator(0x0, rand(Bool,N), rand(Bool,N)) for _ in 1:T]
  phases = 2π*rand(Float64, T)

  println("paulis: ", paulistrings)
  println("phases: ", phases)

  ψ = CAMPS(N)
  arg = (mode == "current") ? collect(1:N) : 0

  for s in 1:3
    P = paulistrings[s]
    ϕ = phases[s]
    C = inv(ψ.Cdag)
    nature = paulinature(arg, C, P)
    println()
    println("=== Step $s P=$P ϕ=$ϕ nature=$nature ===")
    if nature != :trivial
      Q = apply(P, C)
      println("Q (= C * P * C†) = ", Q)
    end
    arg = apply!(ψ, arg, P, ϕ)
    println("After: arg=$arg")
    println("Cdag now: ", ψ.Cdag)

    for j in 1:N
      ex = round(real(expectation(ψ, PauliOperator(0x0, [i==j for i in 1:N], fill(false,N)))), digits=4)
      ey = round(real(expectation(ψ, PauliOperator(0x0, [i==j for i in 1:N], [i==j for i in 1:N]))), digits=4)
      ez = round(real(expectation(ψ, PauliOperator(0x0, fill(false,N), [i==j for i in 1:N]))), digits=4)
      println("  qubit $j: ⟨X⟩=$ex ⟨Y⟩=$ey ⟨Z⟩=$ez")
    end
  end
end

go(ARGS[1])

using QuantumClifford

"After mapping the Pauli string P through the Clifford operator C (ie. P' = C^†PC),
determine the type of action the result has on the state |m⟩^(⊗k) |0⟩^(⊗(N-k)):
- ```:disentanglable```, nontrivial action on at least one |0⟩ qubit:\
action of a phase gate e^iϕP can be incorporated without increasing the MPS bond dimension.
- ```:trivial```, trivial action on all qubits: action of phase gate can be discarded.
- ```:logical```, nontrivial action on |m⟩ qubits only: MPS bond dimension must increase."
function paulinature(free_qubits::Vector{<:Integer}, C::CliffordOperator, P::PauliOperator)
  n = length(P)
  Q = apply(P, C)
  XQ = xbit(Q)
  ZQ = zbit(Q)
  qubits = collect(range(1,n))
  magic_qubits = setdiff(qubits, free_qubits)
  if any(XQ[i] != 0 for i in free_qubits)
    return :disentanglable 
  elseif all(ZQ[i] == 0 for i in magic_qubits) && all(XQ[i] == 0 for i in magic_qubits)
    return :trivial
  else 
    return :logical
  end
end

function paulinature(k::Integer, C::CliffordOperator, P::PauliOperator) 
  N = length(P)
  free = collect(range(k+1, N))
  return paulinature(free, C, P)
end

# === OVERWRITES CliffordMPS ===
# "Append a Clifford gate op acting on the given sites to a Clifford circuit operator C."
# function QuantumClifford.apply!(C::CliffordOperator,
# op::CliffordOperator, sites::AbstractArray{Int,1})
#   apply!(Stabilizer(C.tab), op, sites)
#   return nothing
# end

"Map a Pauli string P through a Clifford operator C, ie. compute C^†PC.\
Warning: not an in-place method!"
function apply(P::PauliOperator, C::CliffordOperator)
  P_tableau = QuantumClifford.Tableau([P])
  apply!(Stabilizer(P_tableau), C)
  return P_tableau[1]
end

"Given a Pauli string P, a Clifford operator C and the free qubits,\
build an analytical disentangling Clifford circuit."
function disentangler(free_qubits::Vector{<:Integer}, C::CliffordOperator, P::PauliOperator)
  Q = apply(P, C)
  i = findfirstfreeXY(Q, free_qubits)
  # println("First free qubit is $i")
  Dtot = one(CliffordOperator, length(Q))

  if zbit(Q)[i]
    Q, phase = reducetoX(Q, i)
    apply!(Dtot, phase)
    # println("Reduced $i")
  end

  Dmain = build_D(Q, i)
  apply!(Dtot, Dmain)

  return Dtot, i
end

function disentangler(k::Integer, C::CliffordOperator, P::PauliOperator)
  N = length(P)
  free = collect(range(k+1, N))
  return disentangler(free, C, P)
end

"Routine for ```disentangler``` to use if the first free qubit with nontrivial action is acted on by a Y"
function reducetoX(P::PauliOperator, i::Integer)
  phasei = one(CliffordOperator, length(P))
  apply!(phasei, tPhase, [i]) # Use 3 phase gates s.t. Y -> X instead of -X
  apply!(phasei, tPhase, [i])
  apply!(phasei, tPhase, [i])
  return apply(P, phasei), phasei
end

findfirstfreeXY(P::PauliOperator, free_qubits::Vector{<:Integer}) = free_qubits[findfirst(i -> (xbit(P)[i] == 1), free_qubits)]

"Routine for ```disentangler``` to use to build the disentangling circuit from Liu and Clark (2025)."
function build_D(Q::PauliOperator, i::Integer)
  tCX = tCNOT
  tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
  tCZ = (tId1 ⊗ tHadamard) * tCNOT * inv(tId1 ⊗ tHadamard)

  D = one(CliffordOperator, length(Q))

  for j in eachindex(Q)
    if j != i
      if Q[j] == (true, false) # X
        apply!(D, tCX, [i, j])
        # println("CX on $i, $j")
      elseif Q[j] == (true, true) # Y
        apply!(D, tCY, [i, j])
        # println("CY on $i, $j")
      elseif Q[j] == (false, true) # Z
        apply!(D, tCZ, [i, j])
        # println("CZ on $i, $j")
      end
    end
  end
  return D
end
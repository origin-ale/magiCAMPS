using QuantumClifford

function paulinature(k::Integer, C::CliffordOperator, P::PauliOperator)
  n = length(P)
  Q = apply(P, C)
  XQ = xbit(Q)
  ZQ = zbit(Q)
  if XQ[k+1:n] != zeros(n-k)
    return :disentanglable # Disentanglable acc. to Liu 2025 thm 1
  elseif ZQ[1:k] == zeros(k) && XQ[1:k] == zeros(k)
    return :trivial
  else 
    return :logical
  end
end

function QuantumClifford.apply!(C::CliffordOperator,
op::CliffordOperator, sites::AbstractArray{Int,1})
  apply!(Stabilizer(C.tab), op, sites)
  return nothing
end

function apply(P::PauliOperator, C::CliffordOperator)
  P_tableau = QuantumClifford.Tableau([P])
  apply!(Stabilizer(P_tableau), C)
  return P_tableau[1]
end

function disentangler(k::Integer, C::CliffordOperator, P::PauliOperator)
  Q = apply(P, C)
  i = findfirstfreeXY(Q, k)
  Dtot = one(CliffordOperator, length(Q))

  if i != k+1
    Q, swap = swapqubits(Q, i, k+1)
    apply!(Dtot, swap)
    i = k+1
  end
  if zbit(Q)[i]
    Q, phase = reducetoX(Q, i)
    apply!(Dtot, phase)
  end
  Dmain = build_D(Q, i)
  apply!(Dtot, Dmain)

  return Dtot
end

function swapqubits(P::PauliOperator, i::Integer, j::Integer)
  swapij = one(CliffordOperator, length(P))
  apply!(swapij, tSWAP, [i,j])
  return apply(P, swapij), swapij
end

function reducetoX(P::PauliOperator, i::Integer)
  phasei = one(CliffordOperator, length(P))
  apply!(phasei, tPhase, [i]) # Use 3 phase gates s.t. Y -> X instead of -X
  apply!(phasei, tPhase, [i])
  apply!(phasei, tPhase, [i])
  return apply(P, phasei), phasei
end

findfirstfreeXY(P::PauliOperator, k::Integer) = findfirst(xbit(P)[k+1:end]) + k

function build_D(Q::PauliOperator, i::Integer)
  tCX = tCNOT
  tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
  tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

  D = one(CliffordOperator, length(Q))

  for j in eachindex(Q)
    if j != i
      if Q[j] == (true, false) # X
        apply!(D, tCX, [i, j])
      elseif Q[j] == (true, true) # Y
        apply!(D, tCY, [i, j])
      elseif Q[j] == (false, true) # Z
        apply!(D, tCZ, [i, j])
      end
    end
  end
  return D
end
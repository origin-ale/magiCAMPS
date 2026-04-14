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
  if i != k+1
    Q = swapqubits!(C, Q, i, k+1) # C modified as a side effect TODO: clean up
    i = k+1
  end
  if zbit(Q)[i]
    Q = reducetoX!(C, Q, i) # C modified as a side effect TODO: clean up
  end
  return build_Ddag(Q, i)
end

function swapqubits!(C::CliffordOperator, P::PauliOperator, i::Integer, j::Integer)
  apply!(C, tSWAP, [i,j])
  swapij = one(CliffordOperator, length(P))
  apply!(swapij, tSWAP, [i,j])
  return apply(P, swapij)
end

function reducetoX!(C::CliffordOperator, P::PauliOperator, i::Integer)
  apply!(C, tPhase*tPhase*tPhase, [i]) # Use 3 phase gates s.t. Y -> X instead of -X
  phasei = one(CliffordOperator, length(P))
  apply!(phasei, tPhase, [i])
  apply!(phasei, tPhase, [i])
  apply!(phasei, tPhase, [i])
  return apply(P, phasei)
end

findfirstfreeXY(P::PauliOperator, k::Integer) = findfirst(xbit(P)[k+1:end]) + k

function build_Ddag(Q::PauliOperator, i::Integer)
  tCX = tCNOT
  tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
  tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

  Ddag = one(CliffordOperator, length(Q))

  for j in eachindex(Q)
    if j != i
      if Q[j] == (true, false) # X
        apply!(Ddag, tCX, [i, j])
      elseif Q[j] == (true, true) # Y
        apply!(Ddag, tCY, [i, j])
      elseif Q[j] == (false, true) # Z
        apply!(Ddag, tCZ, [i, j])
      end
    end
  end
  Ddag = inv(Ddag)
  return Ddag
end
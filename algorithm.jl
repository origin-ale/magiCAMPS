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
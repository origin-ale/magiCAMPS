using ITensors
using ITensorMPS
using QuantumClifford

struct CAMPS # TODO: Change to mutable if immutability becomes annoying
  mps::MPS
  C::Vector
end

"Apply a gate to a quantum circuit with a state encoded by a CAMPS"
function applygate(ψ::CAMPS, gate)
  if isclifford(gate)
    ψp = applyclifford(ψ, gate) 
  else
    if islogical(gate, ψ)
      ψp = applympo(ψ, gate)
    elseif anticommgenerator(gate, ψ.C)
      ψp = addmagic(ψ, gate) # Update both the Clifford and the MPS
    end
  end
  return ψp
end

"""Determine whether a gate is Clifford. \
Currently the program only supports Hadamard, phase, CNOT and T gates."""
function isclifford(gate)
  return QuantumClifford.isclifford(gate)
end

"Given a CAMPS and a Clifford gate, update the CAMPS to one including the new Clifford."
function applyclifford(ψ::CAMPS, gate)
  clifford = copy(ψ.C)
  push!(clifford, gate)
  return CAMPS(ψ.mps, clifford)
end

"Given a Pauli string P and a Clifford gate C, return Pauli string C† P C"
function cliffordconjugate(pauli, gate)
  pauli_tab = Stabilizer([pauli])

  if typeof(gate) <: AbstractTwoQubitOperator
    q1 = gate.q1
    q2 = gate.q2
    gate_cliff = SparseGate(CliffordOperator(typeof(gate)), [q1,q2])
  else
    q = gate.q
    gate_cliff = SparseGate(CliffordOperator(typeof(gate)), [q])
  end

  apply!(pauli_tab, gate_cliff)
  return pauli_tab[1]
end

"Return the number of |0⟩ product states in the MPS"
function nzeros(mps)
  workmps = copy(mps)
  n = 0
  for i in eachindex(mps)
    orthogonalize!(workmps, i)
    t = workmps[i]
    site = siteind(workmps, i)
    densitymtx = t * dag(prime(t, site))
    densitymtx₀₀ = densitymtx[site => 1, prime(site) => 1] 
    (densitymtx₀₀ ≈ 1) && (n += 1)
  end
  return n
end

"Check if a Pauli is a logical operator of CAMPS Clifford."
function islogical(pauli, state)
  cliff = state.C
  n = length(state.mps)
  k = nzeros(state.mps)
  for gate in reverse(cliff)
    pauli = cliffordconjugate(pauli, gate)
  end 
  return xbit(pauli)[n-k+1:end] == zeros(k)
end

"Make an MPO by commuting a Pauli rotation through the Clifford and apply it to the MPS."
function applympo(ψ, rot)
end

"Check if a gate anticommutes with a generator of the Clifford stabilizer."
function anticommgenerator(gate, cliff)
end

"Turn a |0⟩ to an |m⟩ in the product part of the MPS and update the Clifford with the optimal disentangler (Fux 2025)."
function addmagic(ψ, gate)
end
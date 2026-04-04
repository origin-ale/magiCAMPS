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
    if islogical(gate, ψ.C)
      ψp = applympo(ψ, gate)
    elseif anticommgenerator(gate, ψ.C)
      ψp = worstcaseupdate!(ψ, gate) # Update both the Clifford and the MPS
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

"Check if a gate is a logical operator of a Clifford operator."
function islogical(gate, cliff)
end

"Make an MPO by commuting a Pauli rotation through the CAMPS Clifford and apply it to the MPS."
function applympo(ψ, rot)
end

"Check if a gate anticommutes with a generator of the Clifford stabilizer."
function anticommgenerator(gate, cliff)
end

"Update the CAMPS Clifford and Y-rotate a |0⟩ in the MPS."
function worstcaseupdate(ψ, gate)
end
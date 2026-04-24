using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

"```evolve(ψ, t, paulis, phases; [showprogress::Bool])```

Evolve the CAMPS ψ along the Pauli rotation circuit\
with Pauli strings paulis and the given phases, until layer t.
Return the evolved CAMPS."
function evolve(ψ::CAMPS,
                t::Integer,
                paulistrings::Vector{<:PauliOperator},
                phases::Vector{<:Real};
                showprogress = false)
  k = 0
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    k = apply!(ψ, k, paulistrings[s], phases[s], progressbar)
  end
  return ψ, k
end

"```evolve_deepcliffords(ψ, t, paulis, phases; [showprogress::Bool])```

Evolve the CAMPS ψ along a Pauli rotation-doped deep Clifford circuit, until layer t.\
Layer s has a 2N²-deep random Clifford and a paulis[s]-rotation with angle phases[s].
Return the evolved CAMPS."
function evolve_deepcliffords(ψ::CAMPS,
                              t::Integer,
                              paulistrings::Vector{<:PauliOperator},
                              phases::Vector{<:Real};
                              showprogress = false)
  k = 0
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    applyQOp!(ψ, QOp(:randCliffCircuit))
    k = apply!(ψ, k, paulistrings[s], phases[s], progressbar)
  end
  return ψ, k
end

"```apply!(ψ, k, P, ϕ, progressbar)```

Apply the Pauli operator P to the CAMPS ψ with k free qubits, disentangling if possible\
and advancing the progressbar.
Modify ψ in-place and return the new number of free qubits."
function CliffordMPS.apply!(ψ::CAMPS, 
                            k::Integer, 
                            P::PauliOperator, 
                            ϕ::Real, 
                            progressbar::Progress)
  N = length(ψ)
  I = PauliOperator(0x0, fill(false,N), fill(false,N))
  C = inv(ψ.Cdag)

  nature = paulinature(k, C, P)
  if nature == :disentanglable
    ψ.Cdag *= inv(disentangler(k, C, P))
    addmagicstate!(ψ, k, ϕ)
    k += 1
  elseif nature == :logical
    R = PauliSum([cos(ϕ), sin(ϕ)], Stabilizer([I,P]))
    applyGate!(ψ, R)
  elseif nature == :trivial
  end
  next!(progressbar)
  return k
end

"```addmagicstate(ψ, k, phase)```

Turn ψ's (k+1)th qubit from |0⟩ to the Liu and Clark (2025) magic state with given phase"
function addmagicstate!(ψ::CAMPS, k::Integer, phase::Real)
  magifier_os = OpSum()
  magifier_os += cos(phase), "Id", k+1
  magifier_os += -im * sin(phase), "X", k+1
  sites = siteinds(ψ.mps)
  magifier = MPO(magifier_os, sites)
  ψ.mps = ITensors.apply(magifier, ψ.mps)
  return nothing
end
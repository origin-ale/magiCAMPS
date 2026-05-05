using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

bonddim(ψ::CAMPS) = maximum(dim.(linkinds(ψ.mps)))


generate_showvalues(χ, bd) = () -> [("Bond dimension (max $χ)", bd)]

"```evolve_bonddim(ψ, χ, paulis, phases; [showprogress::Bool])```

Evolve the CAMPS ψ along the Pauli rotation circuit\
with Pauli strings paulis and the given phases, until end or bond dim = χ.
Return the evolved CAMPS and stopping time."
function evolve_bonddim(ψ::CAMPS,
                χ::Integer,
                paulistrings::Vector{<:PauliOperator},
                phases::Vector{<:Real};
                showprogress = false)
  free_qubits = collect(range(1,length(ψ)))
  s = 0
  progressthresh = ProgressUnknown(0; dt = 0.05, desc = "Evolving CAMPS… t =", enabled = showprogress)
  while s < length(paulistrings) && bonddim(ψ) < χ
    s += 1
    free_qubits = apply!(ψ, free_qubits, paulistrings[s], phases[s])
    next!(progressthresh; showvalues = generate_showvalues(χ, bonddim(ψ)))
  end
  finish!(progressthresh)
  return ψ, free_qubits, s
end

"```evolve(ψ, t, paulis, phases; [showprogress::Bool])```

Evolve the CAMPS ψ along the Pauli rotation circuit\
with Pauli strings paulis and the given phases, until layer t.
Return the evolved CAMPS."
function evolve(ψ::CAMPS,
                t::Integer,
                paulistrings::Vector{<:PauliOperator},
                phases::Vector{<:Real};
                showprogress = false)
  free_qubits = collect(range(1,length(ψ)))
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    free_qubits = apply!(ψ, free_qubits, paulistrings[s], phases[s])
    next!(progressbar)
  end
  return ψ, free_qubits
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
    k = apply!(ψ, k, paulistrings[s], phases[s])
    next!(progressbar)
  end
  return ψ, k
end

"```apply!(ψ, k, P, ϕ)```

Apply the Pauli operator P to the CAMPS ψ with k free qubits, disentangling if possible.
Modify ψ in-place and return the new number of free qubits."
function CliffordMPS.apply!(ψ::CAMPS, 
                            free_qubits::Vector{<:Integer}, 
                            P::PauliOperator, 
                            ϕ::Real)
  N = length(ψ)
  I = PauliOperator(0x0, fill(false,N), fill(false,N))
  C = inv(ψ.Cdag)

  nature = paulinature(free_qubits, C, P)
  if nature == :disentanglable
    new_Cdag, i = disentangler(free_qubits, C, P)
    ψ.Cdag *= inv(new_Cdag)
    addmagicstate!(ψ, i, ϕ/2) 
    # println("Disentangled")
    deleteat!(free_qubits, findfirst(x -> x == i, free_qubits))
  elseif nature == :logical
    R = PauliSum([cos(ϕ/2), sin(ϕ/2)], Stabilizer([I,P]))
    applyGate!(ψ, R)
    # println("Applied as MPO")
  elseif nature == :trivial
    # println("Trivial")
  end
  # C = inv(ψ.Cdag)
  # @show C
  return free_qubits
end

"```addmagicstate!(ψ, k, phase)```

Turn ψ's given qubit from |0⟩ to the Liu and Clark (2025) magic state with given phase"
function addmagicstate!(ψ::CAMPS, i::Integer, phase::Real)
  magifier_os = OpSum()
  magifier_os += cos(phase), "Id", i
  magifier_os += -im * sin(phase), "X", i
  sites = siteinds(ψ.mps)
  magifier = MPO(magifier_os, sites)
  ψ.mps = ITensors.apply(magifier, ψ.mps)
  return nothing
end
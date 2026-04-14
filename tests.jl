include("algorithm.jl")
using Test

@testset "apply!(clifford, op)" begin
  C0 = one(CliffordOperator, 3) # Identity

  C = deepcopy(C0)
  @assert C == C"XII IXI IIX ZII IZI IIZ" "Initialization or Clifford macro not working as expected!"
  apply!(C, tHadamard, [1])
  @test C == C"ZII IXI IIX XII IZI IIZ"
  apply!(C, tCNOT, [1,3])
  @test C == C"ZII IXI IIX XIX IZI ZIZ"
  apply!(C, tCNOT, [2,1])
  @test C == C"ZZI XXI IIX XIX IZI ZZZ"
  apply!(C, tPhase, [2])
  @test C == C"ZZI XYI IIX XIX IZI ZZZ"
  apply!(C, tHadamard, [2])
  @test C == C"ZXI -XYI IIX XIX IXI ZXZ"

  C = deepcopy(C0)
  @test apply!(C, tHadamard, [2]) === nothing
  @test C != C0

  C = deepcopy(C0)
  apply!(C, tHadamard, [2])
  apply!(C, tHadamard, [2])
  @test C == C0

  C = deepcopy(C0)
  apply!(C, tCNOT, [1, 3])
  apply!(C, tCNOT, [1, 3])
  @test C == C0

  C = deepcopy(C0)
  for _ in 1:4
    apply!(C, tPhase, [2])
  end
  @test C == C0

  op = one(CliffordOperator, 2)
  apply!(op, tHadamard, [1])
  apply!(op, tCNOT, [1, 2])

  C1 = deepcopy(C0)
  apply!(C1, op, [3, 1])

  C2 = deepcopy(C0)
  apply!(C2, tHadamard, [3])
  apply!(C2, tCNOT, [3, 1])

  @test C1 == C2

  @test_throws BoundsError apply!(deepcopy(C0), tHadamard, [0])
  @test_throws BoundsError apply!(deepcopy(C0), tHadamard, [1, 2])

  # Current behavior in QuantumClifford does not reject these malformed indices.
  @test_broken try
    apply!(deepcopy(C0), tHadamard, [4])
    false
  catch e
    e isa BoundsError
  end

  @test_broken try
    apply!(deepcopy(C0), tCNOT, [1])
    false
  catch e
    e isa BoundsError
  end
end

@testset "apply(pauli, cliff)" begin
  C = C"-XX iZY -XZ IZ"
  @test apply(P"YY", C) == P"+ZZ"

  Cid = one(CliffordOperator, 3)
  @test apply(P"YZI", Cid) == P"+YZI"

  C1 = one(CliffordOperator, 1)
  apply!(C1, tPhase, [1])
  @test apply(P"X", C1) == P"+Y"
  @test apply(P"Y", C1) == P"-X"
  @test apply(P"Z", C1) == P"+Z"

  C2 = one(CliffordOperator, 2)
  apply!(C2, tHadamard, [1])
  @test apply(P"XZ", C2) == P"+ZZ"
  @test apply(P"ZZ", C2) == P"+XZ"

  C3 = one(CliffordOperator, 2)
  apply!(C3, tCNOT, [1, 2])
  @test apply(P"XI", C3) == P"+XX"
  @test apply(P"IZ", C3) == P"+ZZ"
  @test apply(P"YI", C3) == P"+YX"
  @test apply(P"IY", C3) == P"+ZY"

  C4 = one(CliffordOperator, 3)
  apply!(C4, tHadamard, [3])
  apply!(C4, tCNOT, [3, 1])
  apply!(C4, tPhase, [1])
  @test apply(P"XIX", C4) == P"+YIZ"
  @test apply(P"ZZI", C4) == P"+ZZZ"
  @test apply(P"IYZ", C4) == P"+YYX"

  @test_throws DimensionMismatch apply(P"XYZ", one(CliffordOperator, 2))
end

@testset "paulinature on C=I" begin
  dnt = :disentanglable
  log = :logical
  trv = :trivial

  C = one(CliffordOperator, 3)
  @test paulinature(0, C, P"IXI") == dnt
  @test paulinature(1, C, P"IXI") == dnt
  @test paulinature(2, C, P"IXI") == log

  @test paulinature(0, C, P"IYI") == dnt
  @test paulinature(1, C, P"IYI") == dnt
  @test paulinature(2, C, P"IYI") == log

  @test paulinature(0, C, P"IZI") == trv
  @test paulinature(1, C, P"IZI") == trv
  @test paulinature(2, C, P"IZI") == log

  @test paulinature(0, C, P"XXI") == dnt
  @test paulinature(1, C, P"XXI") == dnt
  @test paulinature(2, C, P"XXI") == log

  @test paulinature(0, C, P"XZI") == dnt
  @test paulinature(1, C, P"XZI") == log
  @test paulinature(2, C, P"XZI") == log

  @test paulinature(0, C, P"XIZ") == dnt
  @test paulinature(1, C, P"XIZ") == log
  @test paulinature(2, C, P"XIZ") == log
end

@testset "paulinature on C≠I" begin
  dnt = :disentanglable
  log = :logical
  trv = :trivial

  # C1 = H on qubit 1: XZI -> ZZI and ZXI -> XXI.
  C1 = one(CliffordOperator, 3)
  apply!(C1, tHadamard, [1])

  @test paulinature(0, C1, P"XZI") == trv
  @test paulinature(1, C1, P"XZI") == log
  @test paulinature(2, C1, P"XZI") == log

  @test paulinature(0, C1, P"ZXI") == dnt
  @test paulinature(1, C1, P"ZXI") == dnt
  @test paulinature(2, C1, P"ZXI") == log

  # C2 = S on qubits 1 and 3: XIX -> YIY and ZZZ -> ZZZ.
  C2 = one(CliffordOperator, 3)
  apply!(C2, tPhase, [1])
  apply!(C2, tPhase, [3])

  @test paulinature(0, C2, P"XIX") == dnt
  @test paulinature(1, C2, P"XIX") == dnt
  @test paulinature(2, C2, P"XIX") == dnt

  @test paulinature(0, C2, P"ZZZ") == trv
  @test paulinature(1, C2, P"ZZZ") == log
  @test paulinature(2, C2, P"ZZZ") == log

  # C3 = CNOT 1->2: XIZ -> XXZ and IZI -> ZZI.
  C3 = one(CliffordOperator, 3)
  apply!(C3, tCNOT, [1, 2])

  @test paulinature(0, C3, P"XIZ") == dnt
  @test paulinature(1, C3, P"XIZ") == dnt
  @test paulinature(2, C3, P"XIZ") == log

  @test paulinature(0, C3, P"IZI") == trv
  @test paulinature(1, C3, P"IZI") == log
  @test paulinature(2, C3, P"IZI") == log

  # Edge case k = n: upper region is empty.
  Cid = one(CliffordOperator, 3)
  @test paulinature(3, Cid, P"III") == trv
  @test paulinature(3, Cid, P"ZZZ") == log
end

@testset "swapqubits!" begin
  C = one(CliffordOperator, 3)

  @test swapqubits!(C, P"XYZ", 1, 3) == P"ZYX"
  @test C == C"IIX IXI XII IIZ IZI ZII"

  @test swapqubits!(C, P"IIX", 2, 3) == P"IXI"
  @test C == C"IXI IIX XII IZI IIZ ZII"
end

@testset "reducetoX!" begin
  C = one(CliffordOperator, 3)

  @test reducetoX!(C, P"XYZ", 2) == P"XXZ"
  @test C == C"XII -IYI IIX ZII IZI IIZ"

  @test reducetoX!(C, P"YYY", 1) == P"XYY"
  @test C == C"-YII -IYI IIX ZII IZI IIZ"
end
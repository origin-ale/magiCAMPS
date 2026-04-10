include("algorithm.jl")

using Test

N = 6

@testset "applyclifford(ψ, $(split(string(gate))[1]))" for gate in (sHadamard(1), sPhase(2), sCNOT(2,4))
  sites = siteinds("Qubit", N)
  startmps = MPS(sites, "0")
  startcliff = [sHadamard(1), sCNOT(2,3), sPhase(5), sCNOT(1, 6), sHadamard(5)]
  start = CAMPS(startmps, startcliff)
  @test begin
    ψ = deepcopy(start)
    directapply_cliff = [ψ.C..., gate]
    applyclifford(ψ, gate).C == directapply_cliff
  end
end

@testset "isclifford" begin
  @test isclifford(sHadamard(1))
  @test isclifford(sCNOT(2,5))
  @test isclifford(sPhase(6))
  @test !isclifford(sT(5))
end

@testset "cliffordconjugate" begin
  @test cliffordconjugate(P"III", sHadamard(1)) == P"III"
  @test cliffordconjugate(P"III", sCNOT(1,3)) == P"III"
  @test cliffordconjugate(P"III", sPhase(2)) == P"III"
  
  @test cliffordconjugate(P"XYZ", sHadamard(1)) == P"ZYZ"
  @test cliffordconjugate(P"XYI", sCNOT(1,3)) == P"XYX"
  @test cliffordconjugate(P"ZXZ", sPhase(2)) == P"ZYZ"
end

@testset "nzeros" begin
  sites = siteinds("Qubit", N)
  startmps = MPS(sites, "0")
  ψ = deepcopy(startmps)
  @test nzeros(ψ) == N

  ψ = apply(op("H", sites[1]), ψ)
  @test nzeros(ψ) == N-1

  ψ = apply(op("H", sites[1]), ψ)
  @test nzeros(ψ) == N

  ψ = apply(op("CNOT", sites[1], sites[3]), ψ)
  @test nzeros(ψ) == N

  ψ = apply(op("H", sites[1]), ψ)
  @test nzeros(ψ) == N-1

  ψ = apply(op("CNOT", sites[1], sites[3]), ψ)
  @test nzeros(ψ) == N-2
  
   ψ = apply(op("X", sites[N]), ψ)
  @test nzeros(ψ) == N-3
end

@testset "islogical" begin
  sites = siteinds("Qubit", N)
  startmps = MPS(sites, "0")
  startmps[1] = state(sites[1], "1")
  startmps[2] = state(sites[1], "+")
  startcliff = [sHadamard(1), sCNOT(1,2), sCNOT(2,6), sPhase(5)]
  start = CAMPS(startmps, startcliff)

  @test !islogical(P"XXXZZZ", start) broken=true

end
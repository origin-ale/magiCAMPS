include("algorithm.jl")

using Test

N = 6
sites = siteinds("Qubit", N)
startmps = MPS(sites, "0")
startcliff = [sHadamard(1), sCNOT(2,3), sPhase(5), sCNOT(1, 6), sHadamard(5)]

start = CAMPS(startmps, startcliff)

@testset "applyclifford(ψ, $(split(string(gate))[1]))" for gate in (sHadamard(1), sPhase(2), sCNOT(2,4))
  @test begin
    ψ = deepcopy(start)
    directapply_cliff = [ψ.C..., gate]
    applyclifford(ψ, gate).C == directapply_cliff
  end
end

@testset "isclifford($(split(string(gate))[1]))" for gate in (sHadamard(1), sPhase(2), sCNOT(5,6))
  @test isclifford(gate)
end

@testset "isclifford($(split(string(gate))[1]))" for gate in (sT(5),)
  @test !isclifford(gate)
end

@testset "cliffordconjugate" begin
  @test cliffordconjugate(P"III", sHadamard(1)) == P"III"
  @test cliffordconjugate(P"III", sCNOT(1,3)) == P"III"
  @test cliffordconjugate(P"III", sPhase(2)) == P"III"
  
  @test cliffordconjugate(P"XYZ", sHadamard(1)) == P"ZYZ"
  @test cliffordconjugate(P"XYI", sCNOT(1,3)) == P"XYX"
  @test cliffordconjugate(P"ZXZ", sPhase(2)) == P"ZYZ"
end
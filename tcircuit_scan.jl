include("algorithm.jl")
include("evolve_tcircuit.jl")

using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter
using Plots
using LaTeXStrings

Ns = length(ARGS) >= 1 ? parse.(Int, split(ARGS[1], ",")) : [12, 16, 24]
n_samples = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 50
fractions = 0.1:0.1:1.1

results_ee = Dict{Int, Vector{Float64}}()
results_sre = Dict{Int, Vector{Float64}}()

generate_showvalues(toN, sample) = () -> [("t/N",toN), ("sample", sample)]

for N in Ns
	xbits = fill(false, N)
	zbits = fill(false, N)
	zbits[1] = true
	global Z₁ = PauliOperator(0x0, xbits, zbits)
	global T₁ = PauliSum(QOp(:T), N, 1)

	ee_vals = Float64[]
	sre_vals = Float64[]

	progressbar = Progress(length(fractions)*n_samples; desc = "N=$N")
	for f in fractions
		t = round(Int, f * N)

		ee_sum = 0.0
		sre_sum = 0.0
		for sample in 1:n_samples
			ψ = CAMPS(N)
			ψ_evo, k = evolve_tcircuit(ψ, t)

			ee_sum += maximum(eEntropys!(ψ_evo.mps)) / N
			sre_sum += sEntropy(ψ_evo.mps, N^2; α =2) / N
			next!(progressbar, showvalues=generate_showvalues(f, sample))
		end
		avg_ee = ee_sum / n_samples
		avg_sre = sre_sum / n_samples

		push!(ee_vals, avg_ee)
		push!(sre_vals, avg_sre)
	end

	results_ee[N] = ee_vals
	results_sre[N] = sre_vals
end

open("results.dat", "w") do io
	for (i, N) in enumerate(Ns)
		println(io, "# N = $N")
		println(io, "# t/N avg_EE avg_SRE")
		for (j, f) in enumerate(fractions)
			println(io, "$f $(results_ee[N][j]) $(results_sre[N][j])")
		end
		i < length(Ns) && print(io, "\n\n")
	end
end

xs = collect(fractions)
colors = reverse(cgrad(:viridis, length(Ns), categorical = true))

p1 = plot(xlabel = L"t/N", ylabel = L"S_E / N",
          title = "Entanglement entropy vs t/N", legend = :topleft, dpi = 300)
for (i, N) in enumerate(Ns)
	plot!(p1, xs, results_ee[N], label = "N=$N", marker = :circle, color = colors[i])
end
savefig(p1, "avg_ee.png")
display(p1)

p2 = plot(xlabel = L"t/N", ylabel = L"\mathcal{M}_2 / N",
          title = "Rényi stabilizer 2-entropy vs t/N", legend = :topleft, dpi = 300)
for (i, N) in enumerate(Ns)
	plot!(p2, xs, results_sre[N], label = "N=$N", marker = :circle, color = colors[i])
end
savefig(p2, "avg_sre.png")
display(p2)

println("Saved plots to avg_ee.png and avg_sre.png, data to results.dat")

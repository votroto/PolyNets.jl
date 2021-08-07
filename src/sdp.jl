using TKMPMeasures
using SumOfSquares
using MosekTools
using JuMP

function hierarchy(u, S; order=maximum(maxdegree.(u)), N=axes(u, 1))
	m = SOSModel(Mosek.Optimizer)

	@variable(m, w[N])
	@variable(m, μ[i in N], TKMPMeasure; on=S[i], maxdegree=order)

	@objective(m, Min, sum(w))

	@constraint(m, [i in N], ∫(u[i], μ[N .!= i]) <= w[i], domain=S[i], maxdegree=order)
	@constraint(m, [i in N], ∫(1, μ[i]) == 1)

	set_silent(m)
	optimize!(m)

	value.(w), value.(μ)
end
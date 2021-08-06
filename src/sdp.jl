using SumOfSquares
using MosekTools
using JuMP

include("measure_variable.jl")

function hierarchy(u, S; order=maximum(maxdegree.(u)), N=axes(u, 1))
	m = SOSModel(Mosek.Optimizer)

	@variable(m, w[N])
	@variable(m, μ[i in N], MeasureVar; vars=variables(S[i]), maxdegree=order)

	@objective(m, Min, sum(w))

	@constraint(m, [i in N], ∫(u[i], μ[N .!= i]) <= w[i], domain=S[i], maxdegree=order)
	@constraint(m, [i in N], μ[i] >= 0, domain = S[i])
	@constraint(m, [i in N], ∫(1, μ[i]) == 1)

	set_silent(m)
	optimize!(m)

	value.(w), extract_atoms.(μ)
end

using DynamicPolynomials
@polyvar x y z
u1 = 2*x*y^2-x^2-y-z
u2 = -u1
S1 = @set x^2 <= 1 && z^2<=1
S2 = @set 1 - y^2 >= 0

hierarchy([u1, u2], [S1, S2], order=2)
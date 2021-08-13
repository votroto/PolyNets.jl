using MultivariatePolynomials
using TKMPMeasures
using SumOfSquares
using JuMP

function sep∫(p, mus)
	# Exploit separability. Otherwise we'd have to do pre-processing to get rid
	# of quadratic terms induced by the 0th moments.
	(⊂)(X, Y) = all(X .∈ Ref(Y))
	_inner(term) = ∫(term, filter(m -> m.vars ⊂ effective_variables(term), mus))

	sum(_inner, terms(p))
end

function hierarchy(u, S; optimizer, order=maximum(maxdegree.(u)), N=axes(S, 1), silent=true)
	m = SOSModel(optimizer)

	@variable(m, w[N])
	@variable(m, μ[i in N], TKMPMeasure; on=S[i], maxdegree=order)

	@objective(m, Min, sum(w))

	@constraint(m, [i in N], sep∫(u[i], μ[N.!=i]) <= w[i], domain=S[i], maxdegree=order)
	@constraint(m, [i in N], ∫(1, μ[i]) == 1)

	set_silent(m)
	optimize!(m)

	value.(w), value.(μ)
end

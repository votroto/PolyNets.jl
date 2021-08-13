using SumOfSquares
using MosekTools
using JuMP

function measure_expectation(p, m)
	psub(c) = subs(p, (m.variables .=> c.center)...) * c.weight
	sum(psub, m.atoms)
end

function poly_max(u, S; order=maxdegree(u))
	m = SOSModel(Mosek.Optimizer)

	@variable(m, w)
	@objective(m, Min, w)

	@constraint(m, w-u>=0, domain=S, maxdegree=order)

	set_silent(m)
	optimize!(m)

	value(w)
end

function best_response_value(u, S, oponents_strategies; order=maxdegree(u))
	Ep = reduce(measure_expectation, oponents_strategies, init=u)
	poly_max(Ep, S; order)
end
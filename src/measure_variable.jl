using JuMP
using LinearAlgebra
using MultivariatePolynomials
using MultivariateMoments

struct MeasureVarBuilder
	vars
	maxdegree::Int
end

mutable struct MeasureVar
	var_map
	vars
	maxdegree::Int
end

struct NonnegativeOnKConstraintBuilder
	m::MeasureVar
	k::AbstractSemialgebraicSet
end

# -----------------------------------------------------------------------------
# utils
# -----------------------------------------------------------------------------

import Base.-
function -(m::MeasureVar, x::Int64)
	if iszero(x)
		m
	else
		error("Function not implemented")
	end
end

function coeffs(p::AbstractPolynomialLike, vars)
	deg = maxdegree(p)
	ms = monomials(vars, 0:deg) |> reverse
	cs = MultivariatePolynomials.coefficient.(Ref(p), ms, Ref(vars))

	cs, ms
end

function integrate(p::AbstractPolynomialLike, d::MeasureVar)
	trimvec(v::AbstractVector, t) = v[begin:min(t, length(v))]

	l = length(d.var_map.data)
	cs, ms = trimvec.(coeffs(p, d.vars), l)
	us = map(m -> d.var_map[m], ms)

	dot(cs, us)
end

integrate(p::Number, d::MeasureVar) = integrate(p * only(monomials(d.vars, 0)), d)
integrate(p::AbstractPolynomialLike, d::AbstractVector{MeasureVar}) = reduce(integrate, d; init=p)
∫ = integrate

function _moment_matrix(m::MeasureVar)
	xd = monomials(m.vars, 0:m.maxdegree ÷ 2) |> reverse
	Md = xd * xd'
	integrate.(Md, Ref(m))
end

function localizing_matrix(m::MeasureVar, g::AbstractPolynomialLike)
	round_up_even(x) = cld(x,2)*2

	tg = round_up_even(maxdegree(g))
	xd = monomials(m.vars, 0:(m.maxdegree ÷ 2 - maxdegree(g) ÷ 2)) |> reverse
	Ld = xd * xd' * g
	integrate.(Ld, Ref(m))
end

function extract_atoms(m::MeasureVar)
	polyeval(p::AbstractPolynomialLike) = value(p)((variables(p) .=> 1)...)

	b = monomials(m.vars, 0:m.maxdegree ÷ 2) |> reverse
	m = polyeval.(_moment_matrix(m))

	extractatoms(MomentMatrix(m, b), 1e-1)
end

# -----------------------------------------------------------------------------
# JuMP Variables
# -----------------------------------------------------------------------------

function JuMP.build_variable(_err::Function, info::VariableInfo, ::Type{MeasureVar}; vars, maxdegree::Int)
	MeasureVarBuilder(vars, maxdegree)
end

function JuMP.add_variable(model::JuMP.Model, b::MeasureVarBuilder, name::String)
	ms = monomials(b.vars, 0:b.maxdegree) |> reverse
	var_map = @variable(model, [ms]; base_name=name)

	MeasureVar(var_map, b.vars, b.maxdegree)
end

# -----------------------------------------------------------------------------
# JuMP Constraints
# -----------------------------------------------------------------------------

function JuMP.build_constraint(_error::Function, f::MeasureVar, set::MOI.GreaterThan; domain)
	return NonnegativeOnKConstraintBuilder(f, domain)
end

function JuMP.add_constraint(m::Model,	b::NonnegativeOnKConstraintBuilder,	name::String)
	mom_mat = @constraint(m, _moment_matrix(b.m) in PSDCone(), base_name=string(name," moment mat"))
	loc_mats = @constraint(m, [i in inequalities(b.k)], localizing_matrix(b.m, i) in PSDCone(), base_name=string(name," localizing mat ineq"))

	[mom_mat, loc_mats]
end

using MosekTools
using DynamicPolynomials
using SemialgebraicSets
using PolyNets

Base.zero(::Type{Polynomial}) = Polynomial(Int[], PolyVar{true}[]) # dirty

""" trivial two player zero-sum "guessing game" example """
function example_0()
	order = 4
	N = 1:2
	@polyvar x[N]
	S = [@set x[player]^2 <= 1 for player in N]

	u = zeros(Polynomial, length(N), length(N))
	u[1,2] = 2*x[1]*x[2]^2 -x[1]^2 -x[2]
	u[2,1] = -u[1,2]
	usum = sum(u, dims=2)

	ws, ms = hierarchy(usum, S; order, optimizer=Mosek.Optimizer)
	ws, ms, usum, S, order
end

""" Instance of example 2 - a game with a complete graph """
function example_1()
	order=4
	N = 1:3
	@polyvar x[N]
	S = [@set x[player]^2 <= 1 for player in N]

	u = zeros(Polynomial, length(N), length(N))
	u[1,2] = -2*x[1]*x[2]^2 + 5*x[1]*x[2] - x[2]
	u[1,3] = -2*x[1]^2 - 4*x[1]*x[3] - 2*x[3]
	u[2,1] = 2*x[1]*x[2]^2 - 2*x[1]^2 - 5*x[1]*x[2] + x[2]
	u[2,3] = -2*x[2]*x[3]^2 - 2*x[2]^2 + 5*x[2]*x[3]
	u[3,1] = 4*x[1]^2 + 4*x[1]*x[3] + 2*x[3]
	u[3,2] = 2*x[2]*x[3]^2 + 2*x[2]^2 - 5*x[2]*x[3]
	usum = sum(u, dims=2)

	ws, ms = hierarchy(usum, S; order, optimizer=Mosek.Optimizer)
	ws, ms, usum, S, order
end

""" security game example - a proper network game """
function example_3()
	# 2 attackers and 2 defenders
	order = 4
	A, D = 1:2, 3:4
	N = [A; D]
	a = [0.5, 0.8]
	d = [0.7, 0.9]

	# targets
	T = Dict(3=>[1,2], 4=>[3,4])
	Tall = 1:4

	# variables for attackers and defenders
	@polyvar x[A, Tall]
	@polyvar y[D, Tall] # extra variables don't hurt

	# simplex strategy sets
	S = [
		@set x[1,1] + x[1,2] + x[1,3] + x[1,4] == a[1] && x[1,1] >= 0 && x[1,2] >= 0 && x[1,3] >= 0 && x[1,4] >= 0;
		@set x[2,1] + x[2,2] + x[2,3] + x[2,4] == a[2] && x[2,1] >= 0 && x[2,2] >= 0 && x[2,3] >= 0 && x[2,4] >= 0;
		@set y[1,1] + y[1,2] == d[1] && y[1,1] >= 0 && y[1,2] >= 0;
		@set y[2,3] + y[2,4] == d[2] && y[2,3] >= 0 && y[2,4] >= 0;
	]

	# efficiency functions
	e1(x,y) = (1-y^2)*x^2;
	e2(x,y) = (1-y)^2*(2x - x^2);
	e3(x,y) = (1-y^2)*(2x - x^2);
	e4(x,y) = (1-y)^2*x^2;
	e = [e1, e2, e3, e4]

	u = [
		a[1] - sum(e[t](x[1,t], y[j-2,t]) for j in D for t in T[j]) -  sum(a) / length(N);
		a[2] - sum(e[t](x[2,t], y[j-2,t]) for j in D for t in T[j]) -  sum(a) / length(N);
		sum(e[t](x[i,t], y[3-2,t]) for i in A for t in T[3]) -  sum(a) / length(N);
		sum(e[t](x[i,t], y[4-2,t]) for i in A for t in T[4]) -  sum(a) / length(N);
	]

	ws, ms = hierarchy(usum, S; order, optimizer=Mosek.Optimizer)
	ws, ms, u, S, order
end
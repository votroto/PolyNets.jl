
function prettyprint(ws, ms, u, S; N=collect(axes(u,1)))
	println("Game ----------------------------------------------------------- ")
	println("\tplayer set: ", collect(N))
	println("\tstrategy sets:")
	println.(Ref("\t\t"), S)
	println("\tutility functions:")
	println.(Ref("\t\t"), u)
	println("\tNash equilibrium is achieved when ")

	for p in N
		println("\tPlayer ", p, " plays ")

		ch = ["$(round.(a.center, digits=2)) with probability $(round(a.weight, digits=2))" for a in ms[p].atoms]
		println.(Ref("\t\t"), ch)
		println()
	end
	println("\tPayoffs to each player: ", round.(ws, digits=2))
	println(" --------------------------------------------------------------- ")
end
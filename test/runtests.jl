include("examples.jl")
include("best_response.jl")
include("prettyprint.jl")

using Test

function verify_br(ws, ms, ps, ss, o; N=axes(ss, 1))
	for i in N
		br = best_response_value(ps[i], ss[i], ms[N .!= i]; order=o)
		@test isapprox(ws[i], br, atol=2e-2)
	end
end

@testset "examples" begin
	ws, ms, ps, ss, o = example_0()
	verify_br(ws, ms, ps, ss, o)
	prettyprint(ws, ms, ps, ss)

	ws, ms, ps, ss, o = example_1()
	verify_br(ws, ms, ps, ss, o)
	prettyprint(ws, ms, ps, ss)

	ws, ms, ps, ss, o = example_2()
	verify_br(ws, ms, ps, ss, o)
	prettyprint(ws, ms, ps, ss)
end
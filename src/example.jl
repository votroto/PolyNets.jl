include("sdp.jl")

using DynamicPolynomials
@polyvar x y z
u1 = 2*x*y^2-x^2-y-z
u2 = -u1
S1 = @set x^2 <= 1 && z^2<=1
S2 = @set 1 - y^2 >= 0

hierarchy([u1, u2], [S1, S2], order=2)
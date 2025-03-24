# use legendreExpansion to generate first five legendre polynomials
include("legendreExpansion.jl")
using .legendreExpansion

# Generate Legendre polynomials up to degree 5
n = 5
x = -1:0.01:1
legendrePolys = lp_gen_eval(n, x)
# Plot the Legendre polynomials

using Plots
#save the plot to a png file   
# better resolution for the plot
gr(size=(800, 600), dpi=300)

plot(legendrePolys[1], label="P_0(x)", title="Legendre Polynomials", xlabel="x", ylabel="P_n(x)")
plot!(legendrePolys[2], label="P_1(x)")
plot!(legendrePolys[3], label="P_2(x)")
plot!(legendrePolys[4], label="P_3(x)")
plot!(legendrePolys[5], label="P_4(x)")
plot!(legendrePolys[6], label="P_5(x)")
savefig("legendre_polynomials.png")
# I am going to test my legendreExpansion code against the code from LegendrePolynomials.jl and 
# FastGaussQuadrature.jl code by calculating the coefficients of the Legendre Expansion of 1/|x1 - x2|

using FastGaussQuadrature
include("auxFunctions.jl")
using .auxFunctions

n = 100
# r1 > r2
function invDist(x, r1, r2)
    return 1/dist(x, r1, r2)
end

r = 1.0
R = 2.0

# theoretical coefficients r2^l / r1 ^ (l+1)
theo_coeff = [r^l / R^(l+1) for l in 0:n]


# using my own code
include("legendreExpansion.jl")
using .legendreExpansion
x = collect(-1:0.0001:1)
f = invDist.(x, R, r)
@time lp = lp_gen_eval(n, x)
own_coeff = leg_coeff_arr(-1:0.0001:1, invDist.(x, R, r), lp)

# using the code from LegendrePolynomials.jl
using LegendrePolynomials
#=
@time lp2 = collectPl.(x, lmax = n)
r_lp2 = [map(p -> p[i], lp2) for i in firstindex(lp2[1]):lastindex(lp2[1])]
r_lp2 = [collect(v) for v in r_lp2]
pac_coeff = leg_coeff_arr(x, invDist.(x, R, r), r_lp2)
=#
#using LinearAlgebra
#x2, w = gausslegendre(n+1)
#@time eff_coeff = [(2.0 * l + 1.0)/2.0 * dot(w, invDist.(x2, R, r) .* Pl.(x2, l)) for l in 0:n]

using DataFrames, CSV

coeff_df = DataFrame(
    l = 0:n,
    theo_coeff = theo_coeff,
    own_coeff = own_coeff,
    #eff_coeff = eff_coeff
)

    # save the coefficients to a csv file
CSV.write("coefficients.csv", coeff_df)
println(coeff_df)

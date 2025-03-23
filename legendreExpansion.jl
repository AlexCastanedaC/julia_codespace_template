module legendreExpansion
include("mathTools.jl")
using .mathTools

export lp_gen_eval, leg_coeff_arr, series_expansion

# Function that generates Legendre polynomials up to nth degree based on existing array from -1 to 1

function lp_gen_eval(n::Int, x::Vector{Float64})
    legPoly = Vector{Vector{Float64}}(undef, n + 1)
    legPoly[1] = ones(Float64, length(x))
    legPoly[2] = x

    for k in 2:n
        c1 = (2.0 * k - 1) / k
        c2 = (k - 1.0) / k
        legPoly[k + 1] = c1 .* x .* legPoly[k] .- c2 .* legPoly[k - 1]
    end

    return legPoly
end

using Romberg
# The following function creates the coefficients from arrays
function leg_coeff_arr(x::AbstractRange, f::Vector, lP::Vector)
    coefficients = Float64[]
    for i in eachindex(lP)
	    integrand = f .* lP[i]
	    integral = romberg(x, integrand)[1]
	    an = (2.0 * (i-1) + 1.0) * integral / 2.0
	    push!(coefficients,an)
    end
    return coefficients
end

function series_expansion(coefficients, eval_lp)
    return sum([coefficients[i] * eval_lp[i] for i in eachindex(coefficients)])
end

# end of legendreExpansion module
end



# This is the main file for the Julia solver of the spherically-inhomogeneous OZ equation

# We will define our system first, starting with the λ parameter
λ = 0.1
# cutoff radius
R = 10 * λ

# We will define the density function
function ρ0(r::Float64)
	N = 5.0
	V = 4.0 * pi * (R^3) / 3.0
	return N / V
end

include("solNOZ.jl")
using .solNOZ

Δr, r_layers0 = r_layers(5, R)

num_per_layer = n_per_layer(ρ0, r_layers0, Δr)

# number of spherical layers
n = length(r_layers0)

x = range(-1.0, 1.0, length = 20001)
r = Array{Float64, 3}(undef, n, n, length(x))

using LinearAlgebra
# Matrix with number of particles per layer
N = diagm(num_per_layer)

include("auxFunctions.jl")
using .auxFunctions

# three-dimensional array of distances between layers varied by x (cosine of the angle)
for i in 1:n
	for j in 1:n
		r[i, j, :] = dist.(x, r_layers0[i], r_layers0[j])
	end
end

# Now we will choose the potential we will use

include("potentials.jl")
using .potentials

A = 1.0
u = similar(r)
u = yukawa.(r, λ, A)

# we need an appropriate bridge function
# Due to its simplicity, we will use the HNC bridge
include("bridges.jl")
using .bridges
bridge = bridge_HNC(r)

# we will now evaluate the Legendre Polynomials using our x values
include("legendreExpansion.jl")
using .legendreExpansion
n_lP = 50 
lP = lp_gen_eval(n_lP, x)
println("LegPol evaluated")
println("Size of lP: ", size(lP[1]))


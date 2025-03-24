module solNOZ

export r_layers, n_per_layer
# This script contains the necessary functions to solve the OZ inhomogeneous equation
# by using the NG method. The idea is to turn the density ρ(r) into N spherical layers, each containing
# ρ(r) = ρ_i for r_i < r < r_{i+1} where ρ_i is a constant density in the layer i.

# Function to calculate the radius of each layer, which takes the following parameters:
# 1. Number of layers n
# 2. cutoff radius r_cut
function r_layers(n::Int64, r_cut::Float64)
	Δr = r_cut / (2*n)
	r_lay = [(2*i-1)*Δr for i in 1:n]
	return Δr, r_lay
end

# function to calculate number of particles per layer, taking as input the following parameters:
# 1. ρ(r) as a function of r (which will be defined in the main.jl file)
# 2. r_layers array
# 3. Δr

using Romberg

function n_per_layer(ρ::Function, r_layers::Vector{Float64}, Δr::Float64)
	n = length(r_layers)
	n_per_layer = zeros(n)

	# for the moment, this is arbitrary
	int_spacing = Δr/100
	
	for i in 1:n
		r = (r_layers[i] - Δr):int_spacing:(r_layers[i] + Δr)
		integrand = ρ.(r) .* r.^2
		n_per_layer[i] = 4 * pi * romberg(r, integrand)[1]
	end
	return n_per_layer
end
#=
using .legendreExpansion
include("mathTools.jl")
using .mathTools

function simpleIter(n::Number, n_lP::Number,x::AbstractRange, N::Matrix{Float64}, gamma_0::AbstractArray, u::AbstractArray, bridge::AbstractArray, lP::AbstractArray)
	c_dir = exp.(-1.0 .* u .+ gamma_0 .+ bridge) .- gamma_0 .- 1.0
    c_l_coeff = [leg_coeff_arr(x, c_dir[i, j, :], lP) for i in 1:n, j in 1:n]
	c_l_matrix = zeros(n, n, n_lP + 1)
	for i in 1:n, j in 1:n
		c_l_matrix[i, j, :] .= c_l_coeff[i, j]
	end	

	gamma_l = similar(c_l_matrix)
	Iden = I(n) 
	
	for l in 0:n_lP
		den_l = 1.0 / (2.0 * l + 1)
		c = c_l_matrix[:,:,l+1]
		gamma_l[:,:,l+1] =(inv(Iden - den_l.* c*N) * c) - c
	end

	gamma_new_2d = [series_expansion(gamma_l[i, j, :], lP) for i in 1:n, j in 1:n]
	gamma_new = similar(c_dir)

	for i in 1:n, j in 1:n
		gamma_new[i,j,:] .= gamma_new_2d[i, j]
	end
	return gamma_new
end

function inner_product_M(A::Array, B::Array,x::AbstractRange)
	shape_Array = size(A)
	n = shape_Array[1]
	sum = 0.0
	for i in 1:n
		for j in 1:n
                        sum += inner_product(A[i,j,:], B[i,j,:],x) 
                end
        end
	return sum
end

# Define the function to compute corrections and Matrix D
function compute_D(d_n::Vector, x::AbstractRange)
	#Calculate the differences d
	n = length(d_n)
	d = [d_n[end] - d_n[end-i] for i in 1:(n-1)]

	# Calculate the matrix D
	D = [inner_product_M(d[i], d[j], x) for i in 1:(n-1), j in 1:(n-1)]

	d_vec = [inner_product_M(d_n[end], d[i], x) for i in 1:(n-1)] 
	#println(D)
	#println(d_vec)

	return D, d_vec
end

function ng_solve(n::Int64,n_lP::Int64,x::Vector{Float64}, N::Matrix{Float64}, u::AbstractArray, bridge::AbstractArray, initial_gamma::AbstractArray, legendrePol)
	gamma = Vector{Array{Float64}}()
	g = Vector{Array{Float64}}()
	push!(gamma, initial_gamma)
	i = 0
	tol_Ng = 1000.0
	l=0
	# Append new arrays to gamma using simpleIter
	while tol_Ng > 1e-6
		for j in 1:5
			push!(gamma,simpleIter(n, n_lP,x, N, gamma[i + j],u, bridge, legendrePol))
		end
		# Extend g with the last three elements of gamma
		append!(g, gamma[end-4:end])
    	# Compute differences d_n
		d_n = [g[i + k] .- gamma[i + k] for k in 1:5]
		# println(d_n[2][1,1,:])
		# Compute D and d_vec
		D, d_vec = compute_D(d_n, x)
		# Solve for c using linear algebra
		c = D \ d_vec
		# Calculate sum_g
		sum_g = zeros(size(g[end]))
		for j in 1:4
			sum_g .+= c[j] .* g[end-j]
		end
		# Compute f_n1
		f_n1 = (1 - sum(c)) .* g[end] .+ sum_g
		# Calculate difference d_n2
		d_n2 = f_n1 .- g[end]
		# Update tolerance
		tol_Ng = sqrt(inner_product_M(d_n2, d_n2, x))
		# Update the last element of gamma
		gamma[end] = f_n1
		# Increment counters
		i += 5
		l += 1
		println("l = $l, Tol = $tol_Ng")
    end
	return gamma[end]

end
=#
# end module
end

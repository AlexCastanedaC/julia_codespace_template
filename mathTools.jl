module mathTools

export simps_arr, inner_product

function simps_arr(y::Vector, x::Vector)
	n = length(y) - 1
	n % 2 == 0 || error(" y length (number of intervals) must be odd")
	length(x)-1 == n || error("x and y length must be equal")
	h = (x[end]-x[1])/n
	s = y[1] + 4 * sum(y[2:2:end-1]) + 2 * sum(y[3:2:end-2]) + y[end]
	return h/3.0 *s
end

function inner_product(f1::Vector, f2::Vector, x::Vector)
	integrand = f1 .* f2
	return simps_arr(integrand, x)
end

end
module auxFunctions

export dist
# r1 > r2
function dist(x, r1, r2)
    xi = r2 / r1
    return r1 * (1 + xi^2 - 2 * xi * x)^(1/2)
end

end
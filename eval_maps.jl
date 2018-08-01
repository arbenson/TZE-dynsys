function _kth_largest_algebraic(M::Array{Float64,2}, k::Int64)
    (d, V) = eig(M)
    V = real.(V)    
    j = sortperm(real.(d), rev=true)[k]
    return V[:, j] * sign(V[2, j])
end
function kth_largest_algebraic(k::Int64)
    ret(M::Array{Float64,2}) = _kth_largest_algebraic(M, k)
    return ret
end
function largest_algebraic()
    ret(M::Array{Float64,2}) = _kth_largest_algebraic(M, 1)
    return ret
end

function _kth_smallest_algebraic(M::Array{Float64,2}, k::Int64)
    (d, V) = eig(M)
    V = real.(V)    
    j = sortperm(real.(d))[k]
    return V[:, j] * sign(V[2, j])
end
function kth_smallest_algebraic(k::Int64)
    ret(M::Array{Float64,2}) = _kth_smallest_algebraic(M, k)
    return ret
end
function smallest_algebraic()
    ret(M::Array{Float64,2}) = _kth_smallest_algebraic(M, 1)
    return ret
end

function _kth_largest_magnitude(M::Array{Float64,2}, k::Int64)
    (d, V) = eig(M)
    V = real.(V)    
    j = sortperm(abs.(d), rev=true)[k]
    return V[:, j] * sign(V[2, j])
end
function kth_largest_magnitude(k::Int64)
    ret(M::Array{Float64,2}) = _kth_largest_magnitude(M, k)
    return ret
end
function largest_magnitude()
    ret(M::Array{Float64,2}) = _kth_largest_magnitude(M, 1)
    return ret
end

function _kth_smallest_magnitude(M::Array{Float64,2}, k::Int64)
    (d, V) = eig(M)
    V = real.(V)
    j = sortperm(abs.(d))[k]
    return V[:, j] * sign(V[2, j])
end
function kth_smallest_magnitude(k::Int64)
    ret(M::Array{Float64,2}) = _kth_smallest_magnitude(M, k)
    return ret
end
function smallest_magnitude()
    ret(M::Array{Float64,2}) = _kth_smallest_magnitude(M, 1)
    return ret
end
;

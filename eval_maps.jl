using LinearAlgebra

function safer_sign(v::Vector{Float64})
    sgn = sign(v[2])
    if sgn == 0.0; return 1.0; end
    return sgn
end

function _kth_largest_algebraic(M::Array{Float64,2}, k::Int64)
    F = eigen(M)
    d, V = F.values, F.vectors
    j = sortperm(real.(d), rev=true)[k]
    v = real.(V[:, j])
    return v * safer_sign(v)
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
    F = eigen(M)
    d, V = F.values, F.vectors
    j = sortperm(real.(d))[k]
    v = real.(V[:, j])
    return v * safer_sign(v)
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
    F = eigen(M)
    d, V = F.values, F.vectors
    j = sortperm(abs.(d), rev=true)[k]
    v = real.(V[:, j])
    return v * safer_sign(v)    
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
    F = eigen(M)
    d, V = F.values, F.vectors
    j = sortperm(abs.(d))[k]
    v = real.(V[:, j])
    return v * safer_sign(v)    
end
function kth_smallest_magnitude(k::Int64)
    ret(M::Array{Float64,2}) = _kth_smallest_magnitude(M, k)
    return ret
end
function smallest_magnitude()
    ret(M::Array{Float64,2}) = _kth_smallest_magnitude(M, 1)
    return ret
end

function _closest_in_angle(M::Array{Float64,2}, x::Vector{Float64})
    V = eigen(M).vectors
    angles = abs.(vec(x' * V))
    j = findmax(angles)[2]
    v = real.(V[:, j])
    return v * safer_sign(v)    
end
function closest_in_angle(x::Vector{Float64})
    ret(M::Array{Float64,2}) = _closest_in_angle(M, x)
    return ret
end
;

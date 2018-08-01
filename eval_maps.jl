function kth_largest_algebraic(M::Array{Float64,2}, k::Int64)
    (d, V) = eig(M)
    j = sortperm(real.(d), rev=true)[k]
    return V[:, k] * sign(V[1, k])
end
largest_algebraic(M::Array{Float64,2}) =
    kth_largest_algebraic(M::Array{Float64,2}, 1)

function kth_largest_magnitude(M::Array{Float64,2}, k::Int64)
    M = collapse(T, x)
    (d,V) = eig(M)
    j = sortperm(abs.(d), rev=true)[k]
    return V[:, k] * sign(V[1, k])
end
largest_magnitude(M::Array{Float64,2}) =
    kth_largest_algebraic(M::Array{Float64,2}, 1)

function kth_smallest_algebraic(M::Array{Float64,2}, k::Int64)
    M = collapse(T, x)
    (d,V) = eig(M)
    j = sortperm(real.(d))[k]
    return V[:, k] * sign(V[1, k])
end
smallest_algebraic(M::Array{Float64,2}) =
    kth_smallest_algebraic(M::Array{Float64,2}, 1)

function kth_smallest_magnitude(M::Array{Float64,2}, k::Int64)
    M = collapse(T, x)
    (d,V) = eig(M)
    j = sortperm(abs.(d), rev=true)[k]
    return V[:, k] * sign(V[1, k])
end
smallest_magnitude(M::Array{Float64,2}) =
    kth_smallest_magnitude(M::Array{Float64,2}, 1)

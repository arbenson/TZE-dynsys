function check_dimensions(T::Array{Float64}, x::Vector{Float64})
    if !all(size(T) .== length(x))
        error("Vector and tensor dimensions do not match")
    end
end

function apply(T::Array{Float64,3}, x::Vector{Float64})
    check_dimensions(T, x)
    n = length(x)
    y = zeros(Float64,n)
    for k in 1:n
        y += T[:,:,k] * x * x[k]
    end
    return y
end

function apply(T::Array{Float64,4}, x::Vector{Float64})
    check_dimensions(T, x)
    n = length(x)
    y = zeros(Float64,n)
    for l in 1:n, k in 1:n
        y += T[:,:,k,l] * x * x[k] * x[l]
    end
    return y
end

function apply(T::Array{Float64,5}, x::Vector{Float64})
    check_dimensions(T, x)
    n = length(x)
    y = zeros(Float64,n)
    for r in 1:n, l in 1:n, k in 1:n
        y += T[:,:,k,l,r] * x * x[k] * x[l] * x[r]
    end
    return y
end

function apply(T::Array{Float64,6}, x::Vector{Float64})
    check_dimensions(T, x)
    n = length(x)
    y = zeros(Float64,n)
    for s in 1:n, r in 1:n, l in 1:n, k in 1:n
        y += T[:,:,k,l,r,s] * x * x[k] * x[l] * x[r] * x[s]
    end
    return y
end

function collapse(T::Array{Float64,3}, x::Vector{Float64})
    check_dimensions(T, x)    
    n = length(x)
    Y = zeros(Float64, n, n)
    for k in 1:n
        Y += T[:,:,k] * x[k]
    end
    return Y
end

function collapse(T::Array{Float64,4}, x::Vector{Float64})
    check_dimensions(T, x)    
    n = length(x)
    Y = zeros(Float64,n, n)
    for l in 1:n, k in 1:n
        Y += T[:,:,k,l] * x[k] * x[l]
    end
    return Y
end

function collapse(T::Array{Float64,5}, x::Vector{Float64})
    check_dimensions(T, x)
    n = length(x)
    Y = zeros(Float64,n, n)
    for r in 1:n, l in 1:n, k in 1:n
        Y += T[:,:,k,l,r] * x[k] * x[l] * x[r]
    end
    return Y
end

function collapse(T::Array{Float64,6}, x::Vector{Float64})
    check_dimensions(T, x)
    n = length(x)
    Y = zeros(Float64,n, n)
    for s in 1:n, r in 1:n, l in 1:n, k in 1:n
        Y += T[:,:,k,l,r,s] * x[k] * x[l] * x[r] * x[s]
    end
    return Y
end
;

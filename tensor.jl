type SparseTensor3
    I1::Vector{Int64}
    I2::Vector{Int64}
    I3::Vector{Int64}
    V::Vector{Float64}
    dimension::Int64
end
SparseTensor3(I1::Vector{Int64}, I2::Vector{Int64}, I3::Vector{Int64}, V::Vector{Float64}) =
    SparseTensor3(I1, I2, I3, V,
                  max(maximum(I1), maximum(I2), maxmimum(I3)))

type SparseTensor4
    I1::Vector{Int64}
    I2::Vector{Int64}
    I3::Vector{Int64}
    I4::Vector{Int64}    
    V::Vector{Float64}
    dimension::Int64
end
SparseTensor4(I1::Vector{Int64}, I2::Vector{Int64}, I3::Vector{Int64}, I4::Vector{Int64}, V::Vector{Float64}) =
    SparseTensor4(I1, I2, I3, I4, V,
                  max(maximum(I1), maximum(I2), maxmimum(I3), maximum(I4)))

type SparseTensor5
    I1::Vector{Int64}
    I2::Vector{Int64}
    I3::Vector{Int64}
    I4::Vector{Int64}
    I5::Vector{Int64}
    V::Vector{Float64}
    dimension::Int64
end
SparseTensor5(I1::Vector{Int64}, I2::Vector{Int64}, I3::Vector{Int64}, I4::Vector{Int64}, I5::Vector{Int64}, V::Vector{Float64}) =
    SparseTensor4(I1, I2, I3, I4, I5, V,
                  max(maximum(I1), maximum(I2), maxmimum(I3), maximum(I4), maximum(I5)))

# Data types
const SparseTensor = Union{SparseTensor3,SparseTensor4,SparseTensor5}
const DenseTensor = Array{Float64}
const Tensor = Union{DenseTensor,SparseTensor}
const SpFltMat = SparseMatrixCSC{Float64,Int64}

function tensor_apply(T::Array{Float64,3}, x::Vector{Float64})
    n = length(x)
    y = zeros(Float64,n)
    for k in 1:n; y += T[:,:,k] * x * x[k]; end
    return y
end

function tensor_collapse(T::Array{Float64,3}, x::Vector{Float64})
    n = length(x)
    Y = zeros(Float64,n,n)
    for k in 1:n; Y += T[:,:,k] * x[k]; end
    return Y
end

function check_vector_length(T::SparseTensor, x::Vector{Float64})
    if length(x) != T.dimension
        throw("Vector dimension does not match tensor dimension")
    end
end

function apply(T::SparseTensor3, x::Vector{Float64})
    check_vector_length(T, x)
    y = zeros(Float64, T.dimension)
    for (v, i1, i2, i3) in zip(T.V, T.I1, T.I2, T.I3)
        y[i1] += v * x[i2] * x[i3]
    end
    return y
end

function apply(T::SparseTensor4, x::Vector{Float64})
    check_vector_length(T, x)
    y = zeros(Float64, T.dimension)
    for (v, i1, i2, i3, i4) in zip(T.V, T.I1, T.I2, T.I3, T.I4)
        y[i1] += v * x[i2] * x[i3] * x[i4]
    end
    return y
end

function apply(T::SparseTensor5, x::Vector{Float64})
    check_vector_length(T, x)
    y = zeros(Float64, T.dimension)
    for (v, i1, i2, i3, i4, i5) in zip(T.V, T.I1, T.I2, T.I3, T.I4, T.I5)
        y[i1] += v * x[i2] * x[i3] * x[i4] * x[i5] 
    end
    return y
end

function collapse(T::SparseTensor3, x::Vector{Float64})
    check_vector_length(T, x)    
    V = [v * x[i3] for (v, i3) in zip(T.V, T.I3)]
    n = T.dimension
    return convert(SpFltMat, sparse(T.I1, T.I2, V, n, n))
end

function collapse(T::SparseTensor4, x::Vector{Float64})
    check_vector_length(T, x)    
    V = [v * x[i3] * x[i4] for (v, i3, i4) in zip(T.V, T.I3, T.I4)]
    n = T.dimension
    return convert(SpFltMat, sparse(T.I1, T.I2, V, n, n))
end

function collapse(T::SparseTensor5, x::Vector{Float64})
    check_vector_length(T, x)    
    V = [v * x[i3] * x[i4] * x[i5] for (v, i3, i4, i5) in zip(T.V, T.I3, T.I4, T.I5)]
    n = T.dimension
    return convert(SpFltMat, sparse(T.I1, T.I2, V, n, n))
end
;

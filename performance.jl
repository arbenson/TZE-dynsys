function test_tensor(dim::Int64, order::Int64)
    if order == 3
        T = zeros(Float64, dim, dim, dim)
        for i = 1:dim, j = 1:dim, k = 1:dim
            T[i,j,k] = (-1)^i/i + (-1)^j/j + (-1)^k/k
        end
        return T
    elseif order == 4
        T = zeros(Float64, dim, dim, dim, dim)
        for i = 1:dim, j = 1:dim, k = 1:dim, l = 1:dim
            T[i,j,k,l] = (-1)^i/i + (-1)^j/j + (-1)^k/k + (-1)^l/l
        end
        return T
    elseif order == 5
        T = zeros(Float64, dim, dim, dim, dim, dim)    
        for i = 1:dim, j = 1:dim, k = 1:dim, l = 1:dim, r = 1:dim
            T[i,j,k,l,r] = (-1)^i/i + (-1)^j/j + (-1)^k/k + (-1)^l/l + (-1)^r/r
        end
    elseif order == 6
        T = zeros(Float64, dim, dim, dim, dim, dim, dim)        
        for i = 1:dim, j = 1:dim, k = 1:dim, l = 1:dim, r = 1:dim, s = 1:dim
            T[i,j,k,l,r,s] = (-1)^i/i + (-1)^j/j + (-1)^k/k + (-1)^l/l + (-1)^r/r + (-1)^s/s
        end
        return T        
    else
        error("Must have 3 <= order <= 6")
    end
end

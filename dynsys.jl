function kth_largest_algebraic(M::Array{Float64,2} k::Int64)
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

forward_euler(f::Vector{Float64}, x_curr::Vector{Float64}, h::Float64) = x_curr + h * f 

# Integrator returns new iterate
function TZEDS(T::Tensor, Λ, Integrator;
               x0::Vector{Float64}=normalize(ones(Float64,T.dimension)),
               tol::Float64=1e-6,
               maxiter::Int64=100)
    # Data recording
    evec_hist = zeros(length(x0), maxiter + 1)
    evec_hist[:, 1] = x0
    eval_hist = zeros(Float64, maxiter + 1)
    eval_hist[1] = x0' * apply(T, x0)

    # Derivative of the dynamical system
    derivative(u::Vector{Float64}) = Λ(collapse(T, u)) - u
    
    for iter = 1:maxiter
        x_curr = evec_hist[:, iter]
        x_next = real.(Integrator(deriv(x_curr), x_curr))
        x_next /= norm(x_next, 1)
        evec_hist[:, iter + 1] = x_next
        y = apply(T, x_next)
        eval_hist[iter + 1] = x_next' * y
        # Break if we are close enough to an eigenvalue
        rats = abs.(y ./ x_next)
        if (maximum(rats) - minimum(rats)) / minimum(rats) < tol
            break
        end
    end
    return (eval_hist[1:(iter + 1)], evec_hist[1:(iter + 1)])
end

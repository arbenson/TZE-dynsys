include("tensor.jl")
include("integrators.jl")
include("eval_maps.jl")

"""
TZE_dynsys
-------------
The main function for computing tensor Z eigenvectors with dynamical systems.

Input parameters:
- T::Array{Float64}: The (dense) tensor
- Λ: The mapping function, which takes as input type Array{Float64,2} of size n x n
     and returns a vector of size n.
- Integrator: The numerical integrator, which takes as input a derivative
              function and the current iterate. The derivative function maps a
              vector of length n to a vector of length n. The integrator
              function must return the next iterate, given the current iterate
              and the derivative.

Optional parameters:
- x0::Vector{Float64}: the starting guess
- tol::Float64: the stopping tolerance (relative eigenvalue guess)
- maxiter::Int64: the maximum number of numerical integration steps

Returns:
A tuple (eval_hist::Vector{Float64}, evec_hist::Array{Float64,2}), where
eval_hist is the history of eigenvalue estimates, and each column of evec_hist
are the corresponding eigenvector estimate.
"""
function TZE_dynsys(T::Array{Float64}, Λ, Integrator;
                    x0::Vector{Float64}=normalize(ones(Float64,T.dimension)),
                    tol::Float64=1e-6,
                    maxiter::Int64=100)
    # Data recording
    evec_hist = zeros(Float64, length(x0), maxiter + 1)
    evec_hist[:, 1] = x0
    eval_hist = zeros(Float64, maxiter + 1)
    eval_hist[1] = x0' * apply(T, x0)

    # Derivative of the dynamical system
    derivative(u::Vector{Float64}) = Λ(collapse(T, u)) - u

    iter = 1
    while iter <= maxiter
        x_curr = evec_hist[:, iter]
        x_next = real.(Integrator(derivative, x_curr))
        evec_hist[:, iter + 1] = x_next
        y = apply(T, x_next)
        eval_hist[iter + 1] = x_next' * y
        # Break if we are close enough to an eigenvalue
        rats = abs.(y ./ x_next)
        iter += 1
        if maximum(rats) - minimum(rats) < tol
            break
        end
    end
    return (eval_hist[1:iter], evec_hist[1:iter])
end
;

include("tensor.jl")
include("integrators.jl")
include("eval_maps.jl")

# Integrator returns new iterate

function TZEDS(T::Array{Float64}, Λ, Integrator;
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
;

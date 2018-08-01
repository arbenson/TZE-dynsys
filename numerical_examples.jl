include("dynsys.jl")

using Combinatorics
using PyPlot

# Tensor in Example 3.6 from Kolda and Mayo. "Shifted power method for computing
# tensor eigenpairs." SIMAX, 2011.
function T_36()
    function symtensor(T::Array{Float64})
        S = zeros(T)
        d = ndims(T)
        for p in permutations(1:d)
            S += permutedims(T,p)
        end
    return S
    end
    
    T = zeros(Float64, 3,3,3)
    T[1,2,3] = -0.1790
    T[2,3,3] = 0.1773 / 2
    T[1,1,2] = 0.0516 / 2
    T[1,1,3] = -0.0954 /2
    
    T[1,2,2] = -0.1958 /2
    T[1,3,3] = -0.2676 /2
    T[2,2,2] = 0.3251 /6
    T[2,2,3] = 0.2513 /2
    T[3,3,3] = 0.0338 /6

    T = symtensor(T)
    
    T[1,1,1] = -0.1281
    T[2,2,2] = 0.3251
    T[3,3,3] = 0.0338
    
    return T
end

function plots_36(eigenvalue::Float64)
    srand(1)
    T = T_36()
    maxiter = 30

    all_quotients = []
    Λ = kth_smallest_algebraic(2)
    FE = forward_euler(0.5)
    for i in 1:100
        x0=randn(Float64, size(T)[1])
        quotients, xhist = TZE_dynsys(T, Λ, FE, x0=x0, tol=1e-16, maxiter=maxiter)
        push!(all_quotients, quotients)
    end

    figure()
    fsz = 16
    xlabel("Iteration", fontsize=fsz)
    ylabel("Rayleigh quotient", fontsize=fsz)
    if eigenvalue == 0.0018
        quotients = all_quotients[1]
        plot(1:(length(quotients)-1), quotients[2:end])
        ylim(0, 0.002)
        title("V5 (λ = 0.0018)", fontsize=fsz)
    elseif eigenvalue == 0.0033
        quotients = all_quotients[4]
        plot(1:(length(quotients)-1), quotients[2:end])
        ylim(-0.01, 0.004)
        title("V5 (λ = 0.0033)", fontsize=fsz)
    elseif eigenvalue == 0.2294
        quotients = all_quotients[5]
        plot(1:(length(quotients)-1), quotients[2:end])
        ylim(0, 0.25)
        title("V5 (λ = 0.2294)", fontsize=fsz)
    else
        error("Unkown eigenvalue")
    end
    ax = gca()
    ax[:tick_params]("both", labelsize=fsz, length=5, width=1.5)
    tight_layout()
    savefig(string("ex36-V5-", "$(eigenvalue)"[3:end], ".eps"))
end

# Tensor in Example 4.11 from Cui, Dai, and Nie. "All real eigenvalues of
# symmetric tensors.", SIMAX, 2014.
function T_411(n::Int64)
    A = zeros(Float64, n, n, n)
    for i in 1:n, j in 1:n, k in 1:n
        A[i, j, k] = (-1)^i / i + (-1)^j / j + (-1)^k / k
    end
    return A
end

function plots_411(eigenvalue::Float64)
    srand(1)
    n = 5
    T = T_411(n)

    niter = 30
    h = 0.5

    #=
    Map(M) = Mx(M, x -> sortperm(abs.(x), rev=true), 1)
    x0=rand(Float64,n)
    quotients, xhist = Z_evec_dynsys_forward_euler(T,h,niter,Map,normalize(x0))
    pyplot(size=(125,125))
    plot(0:(length(quotients)-1), -quotients, legend=false)
    xlabel!("Iteration")
    ylabel!("Rayleigh quotient")
    title!("V1 (λ = 9.9779)")
    ylims!(4,10.5)
    gui()
    savefig("V1-411.pdf")
    =#

    Map(M) = Mx(M, x -> sortperm(abs.(x)), 1)
    x0 = rand(Float64, n)    
    quotients, xhist = Z_evec_dynsys_forward_euler(T,h,niter,Map,normalize(x0))
    pyplot(size=(125,125))
    plot(0:(length(quotients)-1), -quotients, legend=false)
    xlabel!("Iteration")
    ylabel!("Rayleigh quotient")
    title!("V2 (λ = 0.0000)")
    ylims!(-0.5,4)
    gui()
    savefig("V2-411.pdf")

    #=
    Map(M) = Mx(M, x -> sortperm(real.(x), rev=true), 1)
    x0 = rand(Float64, n)   
    quotients, xhist = Z_evec_dynsys_forward_euler(T,h,niter,Map,normalize(x0))
    pyplot(size=(125,125))
    plot(0:(length(quotients)-1), quotients, legend=false)
    xlabel!("Iteration")
    ylabel!("Rayleigh quotient")
    title!("V3 (λ = 4.2876)")
    ylims!(-5, 5)
    gui()
    savefig("V3-411.pdf")
    =#
end
;

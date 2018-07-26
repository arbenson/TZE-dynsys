using Combinatorics
using Plots

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

function Z_evec_dynsys_forward_euler(
        T::Array{Float64,3}, h::Float64, niter::Int64=30,
        L=(M -> ((d,V) = eig(M); j = sortperm(real(d))[1];
                 V[:,j]*sign(V[1,j]))),
        x0::Vector{Float64}=normalize(ones(Float64,size(T,1))))
    function deriv(u::Vector{Float64})
        return L(tensor_collapse(T, u)) - u
    end
    evec_hist = zeros(length(x0), niter + 1)
    evec_hist[:, 1] = x0
    eval_hist = [x0' * tensor_apply(T, x0)]
    for i = 1:niter
        x_prev = evec_hist[:, i]
        x_next = x_prev + h * deriv(x_prev)
        evec_hist[:, i + 1] = x_next
        push!(eval_hist, real(x_next' * tensor_apply(T, x_next)))
    end
    return (evec_hist, eval_hist)
end

function Mx(M::Array{Float64, 2}, f, k)
    (d,V) = eig(M)
    V = real(V)
    j = f(d)[k]
    v = V[:, j]
    return v * sign(v[2])
end

function T_36()
    function symtensor(T::Array{Float64})
        S = zeros(T)
        d = ndims(T)
        for p in permutations(1:d)
            S += permutedims(T,p)
        end
        S
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
    n = 3
    niter = 30
    h = 0.5

    Map(M) = Mx(M, x -> sortperm(real.(x)), 2)
    all_quotients = []
    for i in 1:100
        x0=randn(Float64,n)
        xhist, quotients = Z_evec_dynsys_forward_euler(T,h,niter,Map,normalize(x0))
        push!(all_quotients, quotients)
    end
    for (i, qs) in enumerate(all_quotients)
        qend = qs[end]
        println("$i $qend")
    end

    if eigenvalue == 0.0018
        pyplot(size=(175,175))
        xlabel!("Iteration")
        ylabel!("Rayleigh quotient")
        quotients = all_quotients[1]
        plot(0:(length(quotients)-1), quotients, legend=false)
        title!("V5 (λ = 0.0018)")
        ylims!(0, 0.002)
        gui()
        savefig("V5a-36.pdf")
    end
    if eigenvalue == 0.0033
        quotients = all_quotients[4]
        plot(0:(length(quotients)-1), quotients, legend=false)
        title!("V5 (λ = 0.0033)")
        ylims!(0, 0.004)
        gui()
        savefig("V5b-36.pdf")
    end
    if eigenvalue == 0.2294
        quotients = all_quotients[5]
        pyplot(size=(150,150))
        plot(0:(length(quotients)-1), quotients, legend=false)
        xlabel!("Iteration")
        ylabel!("Rayleigh quotient")
        title!("V5 (λ = 0.2294)")
        ylims!(0, 0.25)
        gui()
        savefig("V5c-36.pdf")
    end
end

function T_411(n::Int64)
    A = zeros(Float64, n, n, n)
    for i = 1:n
        for j = 1:n
            for k = 1:n
                A[i, j, k] = (-1)^i / i + (-1)^j / j + (-1)^k / k
            end
        end
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
    xhist, quotients = Z_evec_dynsys_forward_euler(T,h,niter,Map,normalize(x0))
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
    xhist, quotients = Z_evec_dynsys_forward_euler(T,h,niter,Map,normalize(x0))
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
    xhist, quotients = Z_evec_dynsys_forward_euler(T,h,niter,Map,normalize(x0))
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

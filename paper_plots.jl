include("dynsys.jl")

using Combinatorics
using MAT
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

function example_36(eigenvalue::Float64)
    srand(1)
    T = T_36()
    maxiter = 20

    all_quotients = []
    Λ = kth_smallest_algebraic(2)
    FE = forward_euler(0.5)
    for i in 1:100
        x0=randn(Float64, size(T)[1])
        tol = -1.0  # negative to do maximum number of iterations
        quotients, xhist = TZE_dynsys(T, Λ, FE, x0=x0, tol=tol, maxiter=maxiter)
        push!(all_quotients, quotients)
    end

    close()
    figure()
    fsz = 24
    xlabel("Iteration", fontsize=fsz)
    ylabel("Rayleigh quotient", fontsize=fsz)
    if eigenvalue == 0.0018
        plot(all_quotients[1], marker="o", lw=3)
        ylim(-0.004, 0.0025)
        title("V5 (λ = 0.0018)", fontsize=fsz)
    elseif eigenvalue == 0.0033
        plot(all_quotients[4], marker="o", lw=3)
        ylim(-0.09, 0.01)
        title("V5 (λ = 0.0033)", fontsize=fsz)
    elseif eigenvalue == 0.2294
        plot(all_quotients[5], marker="o", lw=3)
        ylim(-0.1, 0.25)
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

function example_411(eigenvalue::Float64)
    T = T_411(5)
    maxiter = 20
    tol = -1.0  # negative to run all of the iterations
    FE = forward_euler(0.5) 

    close()
    figure()
    fsz = 24
    xlabel("Iteration", fontsize=fsz)
    ylabel("Rayleigh quotient", fontsize=fsz)
    x0 = normalize(ones(Float64, size(T)[1]))
    figname = ""

    if eigenvalue == 9.9779
        quotients, xhist = TZE_dynsys(T, largest_magnitude(), FE, x0=x0, tol=tol, maxiter=maxiter)
        plot(-quotients, marker="o", lw=3)
        title("V1 (λ = 9.9779)", fontsize=fsz)
        ylim(5, 10.5)
        figname = "ex411-V1.eps"
    elseif eigenvalue == 0.0000
        quotients, xhist = TZE_dynsys(T, smallest_magnitude(), FE, x0=x0, tol=tol, maxiter=maxiter)
        plot(quotients, marker="o", lw=3)
        title("V2 (λ = 0.0000)", fontsize=fsz)
        ylim(-5.5, 0.5)
        figname = "ex411-V2.eps"
    elseif eigenvalue == 4.2876
        quotients, xhist = TZE_dynsys(T, largest_algebraic(), FE, x0=x0, tol=tol, maxiter=maxiter)
        plot(quotients, marker="o", lw=3)
        title("V3 (λ = 4.2876)", fontsize=fsz)
        ylim(-6, 5)
        figname = "ex411-V3.eps"
    else
        error("Unkown eigenvalue")
    end
    ax = gca()
    ax[:tick_params]("both", labelsize=fsz, length=5, width=1.5)
    tight_layout()
    savefig(figname)
end

function scalability(order::Int64)
    function get_time(method::String, dim::Int64)
        endstr = "evals-$order-$dim.mat"
        if method == "DS";      return matread("results/DS-$(endstr)")["time"];              end
        if method == "SDP";     return matread("SDP/results/SDP-$(endstr)")["time"];         end
        if method == "SS-HOPM"; return matread("SS-HOPM/results/SS-HOPM-$(endstr)")["time"]; end
    end
    function get_times(method::String, dims::UnitRange{Int64})
        if method == "DS";      return [get_time("DS",     dim) for dim in dims];  end
        if method == "SDP";     return [get_time("SDP",    dim) for dim in dims];  end
        if method == "SS-HOPM"; return [get_time("SS-HOPM", dim) for dim in dims]; end
        error("Unknown method $method")
    end

    ds_dims = 5:15
    sshopm_dims = 5:15
    sdp_dims = 5:15
    if order == 5; sdp_dims = 5:10; end
    
    ds_times     = get_times("DS",      ds_dims)
    sshopm_times = get_times("SS-HOPM", sshopm_dims)
    sdp_times    = get_times("SDP",     sdp_dims)

    close()
    fsz = 24
    semilogy(collect(ds_dims),     ds_times,     lw=1.5, marker="s", label="DS")
    semilogy(collect(sshopm_dims), sshopm_times, lw=1.5, marker="x", label="SS-HOPM")
    semilogy(collect(sdp_dims),    sdp_times,    lw=1.5, marker="o", label="SDP")
    xlabel("Dimension", fontsize=fsz)
    ylabel("Running time (seconds)", fontsize=fsz)
    if order == 3
        legend(fontsize=fsz-4, loc="upper left", frameon=false,
               labelspacing=0.25)
    end
    title("Order-$(order) tensors", fontsize=fsz)
    ax = gca()
    ax[:tick_params]("both", labelsize=fsz-4, length=6, width=1.5)
    ax[:tick_params]("both", which="minor", length=2, width=1)
    tight_layout()
    savefig("scalability-$(order).eps")
end

function unique_evals(order::Int64)
    function get_evals(method::String, dim::Int64)
        endstr = "evals-$order-$dim.mat"
        evals = Float64[]
        if     method == "DS";      evals = matread("results/DS-$(endstr)")["evals"]
        elseif method == "SDP";     evals = matread("SDP/results/SDP-$(endstr)")["evals"]
        elseif method == "SS-HOPM"; evals = matread("SS-HOPM/results/SS-HOPM-$(endstr)")["evals"]
        else   error("Unknown method $method");
        end
        evals = vec(evals)
        if length(evals) == 0; return evals; end
        # With odd-order tensors, we have sign ambiguity
        if order % 2 == 1; evals = abs.(evals); end
        sort!(evals)
        # Get everything within tolerance 1e-4 (tolerance used by all methods)
        tol = 1e-4
        unique_evals = [evals[1]]
        for eval in evals
            if eval > unique_evals[end] + 1e-4
                push!(unique_evals, eval)
            end
        end
        return unique_evals
    end

    @show get_evals("SDP", 7)
    @show get_evals("DS", 7)
    @show get_evals("SS-HOPM", 7)
end
;

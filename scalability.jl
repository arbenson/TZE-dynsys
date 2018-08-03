include("dynsys.jl")
using MAT

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
        return T
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

function run_test(dim::Int64, order::Int64)
    srand(1234)
    T = test_tensor(dim, order)
    num_trials = 50
    tol = 1e-6
    maxiter = 100

    FE = forward_euler(0.5)
    all_evals = Float64[]
    for k in 1:dim
        for i in 1:num_trials
            x0 = normalize(rand(Float64, dim))
            evals, _, conv = TZE_dynsys(T, kth_largest_algebraic(k), FE, x0=x0, tol=tol, maxiter=maxiter)
            if conv; push!(all_evals, evals[end]); end
        end
        for i in 1:num_trials
            x0 = normalize(rand(Float64, dim))
            evals, _, conv = TZE_dynsys(T, kth_largest_magnitude(k), FE, x0=x0, tol=tol, maxiter=maxiter)
            if conv; push!(all_evals, evals[end]); end
        end
    end
    return all_evals
end

function main()
    for order in 3:5
        println("$(order)...")        
        dimensions = Int64[]
        times = Float64[]
        all_evals = Vector{Vector{Float64}}()
        # Warm-up compilation
        evals = run_test(10, order);
        for dimension in 5:15
            println("\t$(dimension)...")
            tic();
            run_test(dimension, order);
            time = toq();
            # Record data
            matwrite("results/DS-evals-$order-$(dimension).mat",
                     Dict("order"     => order,
                          "dimension" => dimension,
                          "time"      => time,
                          "evals"     => evals))            
        end
    end
end

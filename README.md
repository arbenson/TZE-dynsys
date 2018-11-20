# Tensor Z-eigenvectors & dynamical systems

This code and data repository accompanies the paper

- Computing tensor Z-eigenvectors with dynamical systems. Austin R. Benson and David F. Gleich. [*arXiv:1805.00903*](http://arxiv.org/abs/arXiv:1805.00903), 2018.

All of the code is written in Julia 1.0.

For questions, please email Austin at arb@cs.cornell.edu.

### Setup

This code uses Julia 1.0. It relies on having the following packages installed:

```julia
using Pkg
Pkg.add("Combinatorics")
Pkg.add("FileIO")
Pkg.add("JLD2")
Pkg.add("PyPlot")
```



### Example

Here we show the basics of using the code. The main function is `TZE_dynsys() `. We can load the code in Julia with the following.

```julia
include("dynsys.jl")
```

This function needs 3 inputs:

(1) A dense tensor, represented as type `Array{Float64}`.

Here is an example of a 5-dimensional, fourth-order diagonal tensor.

```julia
dim = 5;
T = zeros(Float64, dim, dim, dim, dim);
for i in 1:dim; T[i, i, i, i] = 2 * (dim + 1 - i); end
```

(2) A function that maps a matrix to one of its eigenvectors. 

We include several maps by default, such as the following (see `eval_maps.jl`).

```julia
# evec for largest algebraic eval
Map1 = largest_algebraic();

# evec for second smallest eval in magnitude
Map2 = kth_smallest_magnitude(2); 

# evec closest to second standard basis vector
Map3 = closest_in_angle(Array(Diagonal(ones(dim))[:,2])); 
```

(3) A numerical integrator function `integrator(f, x)`. The function takes as input a derivative function f and the current iterate x. The derivative function maps a vector to a vector. The integrator function must return the next iterate, given the current iterate and access to the derivative function.

We provide some numerical integrators with the code, so you don't have to worry about all of those details.

```julia
Integrator1 = forward_euler(1.0); # explicit forward Euler with step size 1.0
Integrator2 = forward_euler(0.5); # explicit forward Euler with step size 0.5
Integrator3 = RK4(0.75); # fourth-order explicit Runge-Kutta with step size 0.75
```

With our eigenvector map and numerical integrator in hand, we can now compute Z eigenvectors!

```julia
(evals, evecs, converged) = TZE_dynsys(T, Map3, Integrator2)
```

The vector `evals` contains the Rayleigh quotients at each iteration. The columns of `evecs` are the iterates of the numerical integration scheme. The boolean converged indicates whether or not the numerical integration converged to a tenzor Z-eigenvector.

Here's how we could get all of the eigenvalues for this diagonal tensor.

```julia
FE = forward_euler(0.5)
for k in 1:5
	(evals, evecs, conv) = 
    	TZE_dynsys(T, closest_in_angle(Array(Diagonal(ones(dim))[:,k])), FE); 
	println("converged = $conv");
    println("eval = $(evals[end])");
    println("evec = $(evecs[:,end])");    
end
```



#### Reproduce the figures in the paper

```julia
include("paper_plots.jl")

# Figure 2
example_36(0.0018)  # --> ex36-V5-0018.eps
example_36(0.0033)  # --> ex36-V5-0033.eps
example_36(0.2294)  # --> ex36-V5-2294.eps

# Figure 3
example_411(9.9779) # --> ex411-V1.eps
example_411(0.0000) # --> ex411-V2.eps
example_411(4.2876) # --> ex411-V3.eps

# Figure 4
scalability(3)  # --> scalability-3.eps 
scalability(4)  # --> scalability-4.eps
scalability(5)  # --> scalability-5.eps
```


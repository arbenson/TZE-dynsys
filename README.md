# Tensor Z-eigenvectors & dynamical systems

This code and data repository accompanies the paper

- Computing tensor Z-eigenvectors with dynamical systems. Austin R. Benson and David F. Gleich. [*arXiv:1805.00903*](http://arxiv.org/abs/arXiv:1805.00903), 2018.

All of the code is written in Julia.

For questions, please email Austin at arb@cs.cornell.edu.

## Reproduce the figures and tables in the paper

```julia
include("paper_plots.jl")

plots_36(0.0018)
plots_36(0.0033)
plots_36(0.2294)

plots_411(9.9779)
plots_411(0.0000)
plots_411(4.2876)
```


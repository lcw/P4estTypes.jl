# P4estTypes.jl

[P4estTypes.jl](https://github.com/lcw/P4estTypes.jl) provides a
[Julia](http://julialang.org/) type-based interface to the distributed (via
[MPI](http://www.mpi-forum.org/)) forest-of-octrees library
[p4est](https://p4est.org) by leveraging
[P4est.jl](https://github.com/trixi-framework/P4est.jl) (which provides
low-level Julia bindings to the p4est API).

## Contributing

Contributions are encouraged. In particular only a small part of p4est is
currently exposed. If there are additional functions you would like to use,
please open an [issue](https://github.com/lcw/P4estTypes.jl/issues) or [pull
request](https://github.com/lcw/P4estTypes.jl/pulls).

Additional examples and documentation improvements are also welcome.

## Citation

If you use P4estTypes.jl in your work, please consider citing the following
papers on p4est:
```bibtex
@ARTICLE{BursteddeWilcoxGhattas11,
  author = {Carsten Burstedde and Lucas C. Wilcox and Omar Ghattas},
  title = {{\texttt{p4est}}: Scalable Algorithms for Parallel Adaptive Mesh
           Refinement on Forests of Octrees},
  journal = {SIAM Journal on Scientific Computing},
  volume = {33},
  number = {3},
  pages = {1103-1133},
  year = {2011},
  doi = {10.1137/100791634}
}

@ARTICLE{IsaacBursteddeWilcoxEtAl15,
  author = {Tobin Isaac and Carsten Burstedde and Lucas C. Wilcox and Omar Ghattas},
  title = {Recursive algorithms for distributed forests of octrees},
  journal = {SIAM Journal on Scientific Computing},
  volume = {37},
  number = {5},
  pages = {C497-C531},
  year = {2015},
  doi = {10.1137/140970963}
}
```

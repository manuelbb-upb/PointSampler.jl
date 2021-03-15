```@meta
CurrentModule = PointSampler
```

# PointSampler

This package provides methods to generate samples within a hyper-rectangle in ``ℝ^n``.
It originated from a multiobjective trust region solver that uses local 
surrogate models.

At the moment, there is only one method available, a space-filling, iterative design devised 
by Crombecq in his [dissertation](https://biblio.ugent.be/publication/1970716).
The iterator `MonteCarloThDesign` is described in [Thresholded Monte-Carlo Iterator](@ref).

It is not clear whether more methods will be implemented. \
The monte carlo Voronoi design certainly could be helpful.
The SAMO design is a bit domain specific.

Please see [Surrogates.jl](https://github.com/SciML/Surrogates.jl) for other sampling methods.
It is a more mature package.

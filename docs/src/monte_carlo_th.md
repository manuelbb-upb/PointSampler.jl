# Thresholded Monte-Carlo Iterator

Implementation of a global design to generate samples within a hyper-rectangle 
in ``ℝ^n``. 
It was devised by [Crombecq](https://biblio.ugent.be/publication/1970716) and
is called `mc-intersite-proj-th` in the dissertation.

Samples are generated iteratively to maximize both the intersite distance and the projected
distance.
First, a bunch of samples are generated randomly their scores are calculated.
The projected distance must be greater than a certain **th**reshold.
The best sample is selected to be included in the design.

We provide the iterator `MonteCarloThDesign`.

```@docs
MonteCarloThDesign
```

## Quick Usage Example
This samples 10 points from within ``[0,1]×[0,1]``.

```@example 1
using PointSampler
it = MonteCarloThDesign(; dims = 2, n_points = 10, spawn_factor = 200 );
points = collect(it); 
nothing # hide
```

The points should have distinguishable positions:

```@example 1
using Plots
scatter(Tuple.(points))
savefig("scatter_monte_carlo_th.png"); nothing #hide
```

![](./scatter_monte_carlo_th.png)

### Legacy Function

There is a legacy function to accomplish the same:
```@docs
monte_carlo_th
```
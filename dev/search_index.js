var documenterSearchIndex = {"docs":
[{"location":"monte_carlo_th/#Thresholded-Monte-Carlo-Iterator","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"","category":"section"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"Implementation of a global design to generate samples within a hyper-rectangle  in ℝ^n.  It was devised by Crombecq and is called mc-intersite-proj-th in the dissertation.","category":"page"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"Samples are generated iteratively to maximize both the intersite distance and the projected distance. First, a bunch of samples are generated randomly their scores are calculated. The projected distance must be greater than a certain threshold. The best sample is selected to be included in the design.","category":"page"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"We provide the iterator MonteCarloThDesign.","category":"page"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"MonteCarloThDesign","category":"page"},{"location":"monte_carlo_th/#PointSampler.MonteCarloThDesign","page":"Thresholded Monte-Carlo Iterator","title":"PointSampler.MonteCarloThDesign","text":"MonteCarloThDesign(; dims :: Int ) <: PointIterator\n\nReturn an iterator that provides samples within a hyper- rectangle in ℝ^n where n is the dimension provided  to the iterator as dims.\n\nOptional keyword arguments are: \n\nn_points=100*dims    The maximum number of points to generate.\nlb=zeros(dims)     Lower bounds of hyper-rectangle. Must all be finite.\nub=ones(dims)     Upper bounds of hyper-rectangle. Must all be finite.\nseeds = [].     An initial set of points to add samples to. Is copied and internally modified.\nclean_seeds=true     Throw away seeds that violate the box constraints.\nspawn_factor = 100     In each iteration dims*spawn_factor random points are generated.\nmax_rand_points     Upper bound on the number of random points.\np_dist_th_factor = 0.5   The samples must have a projected distance exceeding 2*p_dist_th_factor/N.\n\n\n\n\n\n","category":"type"},{"location":"monte_carlo_th/#Quick-Usage-Example","page":"Thresholded Monte-Carlo Iterator","title":"Quick Usage Example","text":"","category":"section"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"This samples 10 points from within 0101.","category":"page"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"using PointSampler\nit = MonteCarloThDesign(; dims = 2, n_points = 10, spawn_factor = 200 );\npoints = collect(it); \nnothing # hide","category":"page"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"The points should have distinguishable positions:","category":"page"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"using Plots\nscatter(Tuple.(points))\nsavefig(\"scatter_monte_carlo_th.png\"); # hide\nnothing # hide","category":"page"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"(Image: )","category":"page"},{"location":"monte_carlo_th/#Legacy-Function","page":"Thresholded Monte-Carlo Iterator","title":"Legacy Function","text":"","category":"section"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"There is a legacy function to accomplish the same:","category":"page"},{"location":"monte_carlo_th/","page":"Thresholded Monte-Carlo Iterator","title":"Thresholded Monte-Carlo Iterator","text":"monte_carlo_th","category":"page"},{"location":"monte_carlo_th/#PointSampler.monte_carlo_th","page":"Thresholded Monte-Carlo Iterator","title":"PointSampler.monte_carlo_th","text":"monte_carlo_th( n_points = 10, n_dims = 2; seeds = [], spawn_factor = 50, pdist_threshold_tolerance = 0.5 )\n\nReturn an array of length n_points containing real vectors  representing points in space with n_dims dimensions. The points are iteratively chosen from random point  sets to maximize a space-filling criterion as described in\n\n\"Surrogate Modelling of Computer Experiments with Sequential Experimental Design.\", Crombecq, 2011\n\nThe returned point set is constructed starting with the points in seeds.  If seeds is empty (default), then the singleton set containing the zero vector is used.\n\n\n\n\n\nScale the design returned by the unconstrained version of this  function to the box defined by lb and ub.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = PointSampler","category":"page"},{"location":"#PointSampler","page":"Home","title":"PointSampler","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides methods to generate samples within a hyper-rectangle in ℝ^n. It originated from a multiobjective trust region solver that uses local  surrogate models.","category":"page"},{"location":"","page":"Home","title":"Home","text":"At the moment, there is only one method available, a space-filling, iterative design devised  by Crombecq in his dissertation. The iterator MonteCarloThDesign is described in Thresholded Monte-Carlo Iterator.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is not clear whether more methods will be implemented. \nThe monte carlo Voronoi design certainly could be helpful. The SAMO design is a bit domain specific.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please see Surrogates.jl for other sampling methods. It is a more mature package.","category":"page"}]
}

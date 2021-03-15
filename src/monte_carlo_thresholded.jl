# Space filling designs as in 
# Surrogate Modelling of Computer Experiments with Sequential Experimental Design.", Crombecq, 2011

@doc """
    MonteCarloThDesign(; dims :: Int ) <: PointIterator

Return an iterator that provides samples within a hyper-
rectangle in ``ℝ^n`` where ``n`` is the dimension provided 
to the iterator as `dims`.

Optional keyword arguments are: 

* `n_points=100*dims` 
    The maximum number of points to generate.
* `lb=zeros(dims)`  
    Lower bounds of hyper-rectangle. Must all be finite.
* `ub=ones(dims)`  
    Upper bounds of hyper-rectangle. Must all be finite.
* `seeds = []`.  
    An initial set of points to add samples to. Is copied and internally modified.
* `clean_seeds=true`  
    Throw away seeds that violate the box constraints.
* `spawn_factor = 100`  
    In each iteration `dims*spawn_factor` random points are generated.
* `max_rand_points`  
    Upper bound on the number of random points.
* `p_dist_th_factor = 0.5`
    The samples must have a projected distance exceeding `2*p_dist_th_factor/N`.
"""
@with_kw mutable struct MonteCarloThDesign <: PointIterator
    
    dims :: Int # only mandatory argument

    n_points :: Union{Nothing,Int} = 100 * dims;
    
    lb :: RVec = zeros(dims);
    ub :: RVec = ones(dims);

    seeds :: RVecArr = RVec[];
    clean_seeds :: Bool = true;

    spawn_factor :: Int = 100;
    max_rand_points :: Int = typemax(Int)

    p_dist_th_factor :: Real = 1/2;

    @assert dims > 0
    @assert isnothing(n_points) || ( n_points >= 0)
    @assert length(lb) == length(ub)
    @assert all(isfinite.(lb))
    @assert all(isfinite.(ub))
end

function is_valid( des :: MonteCarloThDesign )
    des.n_points > 0 && des.dims > 0 && ! ( isempty(des.lb) || isempty(des.ub) ) # && n_points == length(point_array) == length(idist_scores) == length( pdist_scores )
end

function clean_seeds!( des :: MonteCarloThDesign )
    if des.clean_seeds
        des.seeds = des.seeds[ good_indices( des.seeds, des.lb, des.ub) ]       
    end
    nothing
end

function scale_seeds!( des :: MonteCarloThDesign ) :: Nothing 
    des.seeds = scale_to_unit_square( des.seeds, des.lb, des.ub )
    nothing
end

function unscale_seeds!( des :: MonteCarloThDesign ) :: Nothing 
    des.seeds = unscale_from_unit_square( des.seeds, des.lb, des.ub )
    nothing
end

function Base.iterate( des :: MonteCarloThDesign )
    if is_valid( des )
        clean_seeds!(des)
        scale_seeds!(des)
        # design has no point associated, initalize with zero vector or first seed
        next_point = isempty( des.seeds ) ? zeros( des.dims ) : des.seeds[1]
        point_array = [ next_point, ]
        return ( unscale_from_unit_square(next_point, des.lb, des.ub), point_array )
    else
        return nothing
    end
end

function Base.iterate( 
    des :: MonteCarloThDesign, point_array :: RVecArr  
    ) :: Union{ Nothing, Tuple{ RVec, RVecArr } }

    n_points_so_far = length(point_array)
    if is_valid( des ) && n_points_so_far < length(des) 
        if n_points_so_far < length( des.seeds )
            # as long as we have fewer points than seeds: 
            # add seed to points
            next_point = des.seeds[ n_points_so_far + 1 ]
        else
            d = des.dims;
            N = n_points_so_far;
            spawn_factor = des.spawn_factor;

            # spawn random candidate points
            candidates = [ rand(d) for j = 1 : min(N * spawn_factor, des.max_rand_points) ];

            # compute candidate scores
            intersite_factor = (( N + 1 )^( 1/d ) - 1)/2
            pdist_factor = (N+1)/2;
            th = 2 * des.p_dist_th_factor / N;

            score_fn = function(p)
                return intersite_factor * minimum( distance( p, point_array ) ) + 
                pdist_factor * minimum( projected_distance_thresholded(p, point_array, th) )
            end
            scores = score_fn.(candidates)
            
            # find point that maximizes score (i.e., is good for intersite distance and 
            # projected distance )
            best_index = argmax( scores )
            next_point = candidates[ best_index ]
        end
        push!( point_array, next_point )

        return unscale_from_unit_square(next_point, des.lb, des.ub), point_array
    end
    return nothing
end

function Base.IteratorSize( des :: MonteCarloThDesign ) 
    if isnothing(des.n_points)
        Base.IsInfinite()
    else
        Base.HasLength()
    end
end
Base.length( des:: MonteCarloThDesign ) = isnothing(des.n_points) ? Inf : des.n_points

## Some legacy functions 
@doc """
    monte_carlo_th( n_points = 10, n_dims = 2; seeds = [], spawn_factor = 50, pdist_threshold_tolerance = 0.5 )

Return an array of length `n_points` containing real vectors 
representing points in space with `n_dims` dimensions.
The points are iteratively chosen from random point 
sets to maximize a space-filling criterion as described in

"Surrogate Modelling of Computer Experiments with Sequential Experimental Design.", Crombecq, 2011

The returned point set is constructed starting with the points in `seeds`. 
If `seeds` is empty (default), then the singleton set containing the zero vector is used.
"""
function monte_carlo_th( 
    n_points :: Int64, n_dims :: Int64 ; 
    seeds :: RVecArr = RVec[],  spawn_factor :: Int64 = 100, 
    pdist_threshold_tolerance :: Real = 0.5 
    )

    iterator = MonteCarloThDesign(;
        dims = n_dims, 
        n_points = n_points, 
        seeds = seeds, 
        spawn_factor = spawn_factor,
        p_dist_th_factor = pdist_threshold_tolerance
    )
    
    return collect( iterator )
end

@doc "Scale the design returned by the unconstrained version of this 
function to the box defined by `lb` and `ub`."
function monte_carlo_th( n_points :: Int64, lb = RVec, ub = RVec; 
    seeds :: RVecArr = RVec[], spawn_factor :: Int64 = 100, 
    pdist_threshold_tolerance :: Real = 0.5, clean_seeds ::Bool = true 
    )

    iterator = MonteCarloThDesign(;
        dims = length(lb), 
        n_points = n_points,
        lb = lb, 
        ub = ub, 
        seeds = seeds, 
        clean_seeds = clean_seeds,
        spawn_factor = spawn_factor,
        p_dist_th_factor = pdist_threshold_tolerance
    )
    
    return collect( iterator )
end

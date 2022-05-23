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
struct MonteCarloThDesign{
VT <: AbstractVector{<:AbstractFloat},
ST <: AbstractVector{<:VT}
} <: PointIterator
    
    dims :: Int # only mandatory argument

    n_points :: Union{Nothing,Int}
    
    lb :: VT 
    ub :: VT 
    width :: VT
    
    seeds :: ST

    spawn_factor :: Int
    max_rand_points :: Int

    p_dist_th_factor :: Float64

    was_scaled :: Bool 
end

function MonteCarloThDesign(; 
    dims :: Int, n_points :: Int = 100 * dims, 
    lb :: AbstractVector{<:Real} = zeros(dims), 
    ub :: AbstractVector{<:Real} = ones(dims), 
    seeds :: AbstractVector{<:AbstractVector{<:Real}} = Vector{Float64}[], 
    clean_seeds :: Bool = true, 
    spawn_factor :: Int = 20, 
    max_rand_points :: Int = typemax(Int), 
    p_dist_th_factor :: Real = 0.5,
    needs_scaling :: Bool = !iszero(lb) && !all(isone.(ub)) 
)

    @assert dims > 0 "`dims` must be positive."
    @assert isnothing(n_points) || n_points >= 0 "`n_points` must be non-negative."
    @assert length(lb) == dims "`dims` and length of variable bounds must be the same."
    @assert length(lb) == length(ub) "Lower bounds `lb` and upper bounds `up` must have same length."
    @assert all(isfinite.(lb)) "All lower bounds `lb` must be finite."
    @assert all(isfinite.(ub)) "All upper bounds `ub` must be finite."
    @assert spawn_factor > 0 "`spawn_factor` must be positive."
    @assert max_rand_points > 0 "`max_rand_points` must be positive."

    FT = Base.promote_eltype( lb, ub, eltype(seeds), Float32 )
    _VT = _vec_type( lb, FT )

    _lb = convert( _VT , lb )
    _ub = convert( _VT, ub )

    w = _ub .- _lb

    _seeds = convert( _vec_type( seeds, _VT ), convert.(_VT , seeds ) )
    _seeds = if clean_seeds
        _seeds[ good_indices( _seeds, _lb, _ub ) ]
    else
        _seeds
    end

    _seeds = if needs_scaling
        scale_with_lb_and_width( _seeds, _lb, w)
    else
        _seeds
    end

    return MonteCarloThDesign(dims, n_points, _lb, _ub, w, _seeds,
        spawn_factor, max_rand_points, Float64(p_dist_th_factor), needs_scaling
    )
end

function unscale_result( x, des :: MonteCarloThDesign )
    des.was_scaled && return unscale_with_lb_and_width( x, des.lb, des.width )
    return x
end

# NOTE / TODO 
#= 
Stateful iterators were introduced after I finished this.
The external state `point_array` could easily be stored in `des`.
=#

"Return the first design point (within original variable boundaries) 
and the state -- a list of points in unit hypersquare generated so far."
function Base.iterate( des :: MonteCarloThDesign{VT,ST} ) where {VT, ST}
    if des.n_points > 0
        # design has no point associated, initalize with zero vector or first seed
        next_point = isempty( des.seeds ) ? convert(VT, zeros(des.dims)) : des.seeds[1]
        point_array = [ next_point, ]
        return ( unscale_result(next_point, des), point_array )
    else
        return nothing
    end
end

function _rand_candidates( VT :: Type, d, N) :: Vector{VT}
    T = eltype(VT)
    return [ VT(rand(T,d)) for _ = 1 : N  ]
end

function score_fn(p :: RVec, point_array, d, N, α)
    intersite_factor = (( N + 1 )^( 1/d ) - 1)/2
    pdist_factor = (N+1)/2
    th = α/N #des.p_dist_th_factor / N

    pdist = minimum( projected_distance_vector(p, point_array) )
    idist = minimum(distance_vector(p, point_array))
    if pdist >= th 
        return (intersite_factor + pdist_factor) * idist
    else
        return intersite_factor * idist
    end
end

function score_fn(P :: RVecArr, args... )
    return [ score_fn( p, args...) for p = P ]
end

function Base.iterate( des :: MonteCarloThDesign{VT,ST}, point_array ) where {VT, ST}

    n_points_so_far = length(point_array)
    if n_points_so_far < length(des) 
        if n_points_so_far < length( des.seeds )
            # as long as we have fewer points than seeds: 
            # add seed to points
            next_point = des.seeds[ n_points_so_far + 1 ]
        else
            d = des.dims
            N = n_points_so_far
            spawn_factor = des.spawn_factor

            # spawn random candidate points
            num_candidates = min(N * spawn_factor, des.max_rand_points)
            candidates = _rand_candidates(VT, d, num_candidates)

            # compute candidate scores
            scores = score_fn(candidates, point_array, d, N, des.p_dist_th_factor )
            
            # find point that maximizes score (i.e., is good for intersite distance and 
            # projected distance )
            best_index = argmax( scores )
            next_point = candidates[ best_index ]
        end
        push!( point_array, next_point )

        return unscale_result(next_point, des), point_array
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
Base.length( des:: MonteCarloThDesign ) :: Union{Float64, Int} = isnothing(des.n_points) ? Inf : des.n_points

# so that collect returns the right kind of array
Base.IteratorEltype( :: MonteCarloThDesign ) = Base.HasEltype()
Base.eltype( :: MonteCarloThDesign{VT,ST} ) where {VT,ST} = VT 

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
    n_points :: Integer, n_dims :: Integer;
    seeds :: RVecArr = Vector{Float64}[], spawn_factor :: Integer = 100, 
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
function monte_carlo_th( n_points :: Integer, lb = RVec, ub = RVec; 
    seeds :: RVecArr = Vector{Float64}[], spawn_factor :: Integer = 100, 
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

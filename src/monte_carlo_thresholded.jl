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
@with_kw struct MonteCarloThDesign{
    VT <: AbstractVector{<:AbstractFloat},
    ST <: AbstractVector{<:VT}
    } <: PointIterator
    
    dims :: Int # only mandatory argument

    n_points :: Union{Nothing,Int} = 100 * dims
    
    lb :: VT = zeros(dims)
    ub :: VT = ones(dims)

    seeds :: ST = Vector{ eltype(lb) }[]
    clean_seeds :: Bool = true

    spawn_factor :: Int = 100
    max_rand_points :: Int = typemax(Int)

    p_dist_th_factor :: Float64 = .5

    needs_scaling :: Bool = !(iszero(lb) && all(ub .== 1))

    function MonteCarloThDesign{VT,ST}( dims :: Int, n_points :: Int, lb :: VT, ub :: VT, seeds :: ST, clean_seeds :: Bool, 
        spawn_factor :: Int, max_rand_points :: Int, p_dist_th_factor, needs_scaling :: Bool ) where{ VT, ST }

        @assert dims > 0 "`dims` must be positive."
        @assert isnothing(n_points) || n_points >= 0 "`n_points` must be non-negative."
        @assert length(lb) == dims "`dims` and length of variable bounds must be the same."
        @assert length(lb) == length(ub) "Lower bounds `lb` and upper bounds `up` must have same length."
        @assert all(isfinite.(lb)) "All lower bounds `lb` must be finite."
        @assert all(isfinite.(ub)) "All upper bounds `ub` must be finite."

        copied_seeds = deepcopy(seeds)
        return new{VT,ST}(dims,n_points,lb,ub,copied_seeds,clean_seeds,spawn_factor,max_rand_points,p_dist_th_factor,needs_scaling)
    end

    # constructor for when vectors with non-float entries are provided
    function MonteCarloThDesign(dims, n_points, lb :: RVec, ub :: RVec, 
        seeds :: RVecArr, clean_seeds, spawn_factor, max_rand_points, p_dist_th_factor, needs_scaling )
        T = Base.promote_eltype( lb, ub, eltype(seeds), Float64 )
        
        # we simply go with the array type of `lb` for `new_lb,new_ub,new_seeds`
        # to avoid abstract common types, e.g. promote_type( Vector{Float64}, SVector{3,Float64}) = AbstractVector{Float64}
        new_lb = similar(lb, T)
        new_lb[:] = [convert(T, x) for x in lb]

        new_ub = similar(lb, T)
        new_ub[:] = [convert(T, x) for x in ub]

        new_seeds = similar(seeds, typeof(new_lb))
        for (i,s) in enumerate(new_seeds)
            s[:] = [ convert(T, x) for x in seeds[i] ]
        end

        return MonteCarloThDesign(dims, n_points, new_lb, new_ub, new_seeds, clean_seeds, spawn_factor, max_rand_points, p_dist_th_factor, needs_scaling )
    end
end

function clean_seeds!( des :: MonteCarloThDesign )
    if des.clean_seeds
        cleaned_seeds = copy(des.seeds[good_indices( des.seeds, des.lb, des.ub)])
        empty!(des.seeds)
        append!(des.seeds, cleaned_seeds)
    end
    return nothing
end

function scale_seeds!( des :: MonteCarloThDesign ) :: Nothing 
    if des.needs_scaling
        for seed in des.seeds 
            scale_to_unit_square!( seed, des.lb, des.ub )
        end
    end
    return nothing
end

function unscale_seeds!( des :: MonteCarloThDesign ) :: Nothing 
    if des.needs_scaling
        for seed in des.seeds
            unscale_from_unit_square!( seed, des.lb, des.ub )
        end
    end
    return nothing
end

function unscale_result( x, des :: MonteCarloThDesign )
    if des.needs_scaling
        return unscale_from_unit_square( x, des.lb, des.ub )
    else
        return x
    end
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
        # throw out seeds that are not in variable bounds
        clean_seeds!(des)
        # scale seeds to unit hypersquare [0,1]^dims
        scale_seeds!(des)
            
        # design has no point associated, initalize with zero vector or first seed
        next_point = isempty( des.seeds ) ? VT(zeros(des.dims)) : des.seeds[1]
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
            intersite_factor = (( N + 1 )^( 1/d ) - 1)/2
            pdist_factor = (N+1)/2
            th = 2 * des.p_dist_th_factor / N

            function (score_fn(p :: AbstractVector)::Float64)
                return intersite_factor * minimum( distance( p, point_array ) ) + 
                    pdist_factor * minimum( projected_distance_thresholded(p, point_array, th) )
            end
            scores = score_fn.(candidates) :: Vector{Float64}
            
            # find point that maximizes score (i.e., is good for intersite distance and 
            # projected distance )
            best_index = argmax( scores ) :: Int
            next_point = candidates[ best_index ] :: VT
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
    pdist_threshold_tolerance :: Float64 = 0.5 
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
    pdist_threshold_tolerance :: Float64 = 0.5, clean_seeds ::Bool = true 
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

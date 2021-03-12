
@doc """
Return the minium projected distance between two point vectors `p1` and `p2`,
i.e. the minimum absolute value of the differences of the coordinates of `p1` and `p2`.
"""
function projected_distance( p1 :: RVec, p2 :: RVec )
    minimum( abs.( p1 .- p2) )
end

@doc "Return array of projected distances for point `p1` against every point in `P`."
function projected_distance( p1 :: RVec, P :: RVecArr )
    [ projected_distance(p1, p2) for p2 ∈ P ]
end

@doc "Return array of projected distances for point `p1` against every point in `P`, but every value below `threshold` is set to `0.0`."
function projected_distance_thresholded( p1 :: RVec, P :: RVecArr, threshold :: Real = 0.1 )
    pdist_array = projected_distance( p1, P );
    pdist_array[ pdist_array .< threshold ] .= 0
    return pdist_array
end

@doc "Return Euclidean distance from `p1` to every point in `P`."
function distance( p1 :: RVec, P :: RVecArr )
    [ norm(p1 .- p2, 2) for p2 ∈ P ]
end

"Indices of vectors in `set` that violate the box constraints `lb` or `ub`."
function bad_indices( set :: RVecArr, lb :: RVec, ub = RVec )
    [ any( s .< lb ) || any( s .> ub ) for s ∈ set]
end

"Indices of vectors in `set` that conform to the box constraints `lb` or `ub`."
function good_indices( set :: RVecArr, lb :: RVec, ub = RVec )
    .!bad_indices(set, lb, ub)
end

"Delete vectors from array `seeds` that violate the box constraints `lb` and `ub`."
function discard_bad_seeds!( seeds :: RVecArr, lb :: RVec, ub = RVec )
    bad_seeds = bad_indices( seeds, lb, ub )
    deleteat!(seeds, bad_seeds)
end

# Simple scaling from and to the square [0,1]ⁿ for finite box constraints 

@memoize function _width( lb :: RVec, ub :: RVec ) :: RVec 
    @assert length(ub) == length(lb)
    @assert all( isfinite.(lb) )
    @assert all( isfinite.(ub) )
    return ub .- lb
end

function scale_to_unit_square( p :: RVec, lb :: RVec, ub :: RVec ) :: RVec
    w = _width(lb, ub)
    return ( p .- lb ) ./ w
end

function scale_to_unit_square!( p :: RVec, lb :: RVec, ub :: RVec ) :: Nothing
    w = _width(lb, ub)
    p .-= lb
    p ./= w;
    nothing
end

function scale_to_unit_square( P :: RVecArr, lb :: RVec, ub :: RVec ) :: RVecArr
    [ scale_to_unit_square( p, lb, ub ) for p ∈ P ]
end

function scale_to_unit_square!( P :: RVecArr, lb :: RVec, ub :: RVec ) :: Nothing
    for p ∈ P 
        scale_to_unit_square!( p, lb, ub )
    end
    nothing
end

function unscale_from_unit_square( p :: RVec, lb :: RVec, ub :: RVec ) :: RVec 
    w = _width(lb, ub)
    return lb .+ w .* p
end

function unscale_from_unit_square!( p :: RVec, lb :: RVec, ub :: RVec ) :: Nothing 
    w = _width(lb, ub)
    p .*= w
    p .+= lb
    nothing
end

function unscale_from_unit_square( P :: RVecArr, lb :: RVec, ub :: RVec ) :: RVecArr
    [ unscale_from_unit_square(p, lb, ub) for p ∈ P ]
end

function unscale_from_unit_square!( P :: RVecArr, lb :: RVec, ub :: RVec ) :: Nothing
    for p ∈ P 
        unscale_from_unit_square!(p, lb, ub)
    end
    nothing
end
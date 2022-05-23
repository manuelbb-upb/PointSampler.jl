
@doc """
Return the minium projected distance between two point vectors `p1` and `p2`,
i.e. the minimum absolute value of the differences of the coordinates of `p1` and `p2`.
"""
function projected_distance( p1 :: RVec, p2 :: RVec )
    minimum( abs.( p1 .- p2) )
end

@doc "Return array of projected distances for point `p1` against every point in `P`."
function projected_distance_vector( p1 :: RVec, P :: RVecArr )
    [ projected_distance(p1, p2) for p2 ∈ P ]
end

@doc "Return Euclidean distance from `p1` to every point in `P`."
function distance_vector( p1 :: RVec, P :: RVecArr )
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

# Simple scaling from and to the square [0,1]ⁿ for finite box constraints 

function _width( lb :: RVec, ub :: RVec )
    #= 
    @assert length(ub) == length(lb)
    @assert all( isfinite.(lb) )
    @assert all( isfinite.(ub) )
    =#
    return ub .- lb
end

function scale_with_lb_and_width( p :: RVec, lb, w )
    return ( p .- lb ) ./ w
end

function scale_to_unit_square( p :: RVec, lb :: RVec, ub :: RVec ) 
    w = _width(lb, ub)
    return scale_with_lb_and_width(p, lb, w)
end

function scale_with_lb_and_width!( p :: RVec, lb, w ) :: Nothing
    p .-= lb
    p ./= w;
    return nothing
end

function scale_to_unit_square!( p :: RVec, lb :: RVec, ub :: RVec ) :: Nothing
    w = _width(lb, ub)
    return scale_with_lb_and_width!(p, lb, w)
end

function scale_with_lb_and_width( P :: RVecArr, lb, w )
    return [ scale_with_lb_and_width(p, lb, w ) for p ∈ P ]
end

function scale_to_unit_square( P :: RVecArr, lb :: RVec, ub :: RVec )
    w = _width(lb, ub)
    return scale_with_lb_and_width(P, lb, w )
end

function scale_with_lb_and_width!( P :: RVecArr, lb, w )
    for p in P 
        scale_with_lb_and_width!(p, lb, w ) 
    end
    return nothing
end

function scale_to_unit_square!( P :: RVecArr, lb :: RVec, ub :: RVec ) :: Nothing
    w = _width(lb, ub)
    return scale_with_lb_and_width!( P, lb, w )
end

function unscale_with_lb_and_width( p :: RVec, lb, w )
    return lb .+ w .* p
end

function unscale_from_unit_square( p :: RVec, lb :: RVec, ub :: RVec ) 
    w = _width(lb, ub)
    return unscale_with_lb_and_width(p, lb, w )
end

function unscale_with_lb_and_width!( p :: RVec, lb, w )
    p .*= w 
    p .+= lb 
    return nothing
end

function unscale_from_unit_square!( p :: RVec, lb :: RVec, ub :: RVec ) 
    w = _width(lb, ub)
    return unscale_with_lb_and_width(p, lb, w )
end

function unscale_with_lb_and_width( P :: RVecArr, lb, w )
    [ unscale_with_lb_and_width(p, lb, w) for p ∈ P ]
end

function unscale_from_unit_square( P :: RVecArr, lb :: RVec, ub :: RVec ) 
    w = _width(lb, ub)
    return unscale_with_lb_and_width(P, lb, ub)
end

function unscale_with_lb_and_width!( P :: RVecArr, lb, w )
    for p in P
        unscale_with_lb_and_width(p, lb, w)
    end
    return nothing
end

function unscale_from_unit_square!( P :: RVecArr, lb :: RVec, ub :: RVec ) :: Nothing
    w = _width(lb, ub)
    return unscale_with_lb_and_width!(P, lb, w)
end

"Return the type of a vector with elements of type `F`."
_vec_type( ::Type{<:Vector}, F ) = Vector{F}
"Return the type of a static vector with elements of type `F`."
_vec_type( T::Type{<:StaticVector}, F ) = similar_type(T, F)
_vec_type( x :: T, F ) where T<:AbstractVector = _vec_type(T, F)

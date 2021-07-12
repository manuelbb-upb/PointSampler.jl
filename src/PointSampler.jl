module PointSampler
using LinearAlgebra: norm
using Memoize: @memoize
using Parameters: @with_kw

export monte_carlo_th, MonteCarloThDesign

const RVec = AbstractVector{<:Real}
const RVecArr = AbstractVector{<:RVec};

abstract type PointIterator end

include("utils.jl")

include("monte_carlo_thresholded.jl");

end
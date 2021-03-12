using PointSampler
using Test

@testset "PointSampler.jl" begin
    include("monte_carlo_th.jl");
    include("monte_carlo_th_iterator.jl");
end

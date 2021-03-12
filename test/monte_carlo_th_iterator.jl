using PointSampler
using Test

P1 = collect(MonteCarloThDesign(10,1))
@test [0.0] ∈ P1
@test length(P1) == 10

P1_2 = collect(MonteCarloThDesign(10, [1.0], [5.0]))
@test [1.0] ∈ P1_2
@test min( P1_2... )[end] >= 1.0
@test max( P1_2... )[end] <= 5.0

seeds = [ rand(2) for i = 1:5 ]
P2_1 = collect(MonteCarloThDesign(3, zeros(2), ones(2), seeds ))
@test length( P2_1 ) == 3
@test all( seeds[1:3] .≈ P2_1 )

P2_2 = collect(MonteCarloThDesign(3, 2.0 .* ones(2), 3.0 .* ones(2), seeds ))   # clean_seeds = true (default)
@test !( seeds ⊆ P2_2 )

push!(seeds, [2.5; 2.5])
P2_3 = collect(MonteCarloThDesign(3, 2.0 .* ones(2), 3.0 .* ones(2), seeds ; clean_seeds = true ))
@test [2.5; 2.5] ∈ P2_3

P2_4 = collect(MonteCarloThDesign(3, 2.0 .* ones(2), 3.0 .* ones(2), seeds ; clean_seeds = false ))
@test all( seeds[1:3] .≈ P2_4 )

P2_5 = collect(MonteCarloThDesign(10, zeros(2), ones(2), seeds))
@test length(P2_5) == 10
@test all( seeds[1:5] .≈ P2_5[1:5] )
@test !( seeds[6] ≈ P2_5[6])    # sixt seed should have been filtered out

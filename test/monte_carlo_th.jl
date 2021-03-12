using PointSampler
using Test

# Unit Line Design (1D) with 10 points, including zero
P1 = monte_carlo_th( 10, 1 )

@test [0.0] ∈ P1
@test length(P1) == 10
@test max( P1... )[end] <= 1.0
@test min( P1... )[end] >= 0.0

# Unit Square Design (2D) with 20 points, including zero vector
P2 = monte_carlo_th( 20, 2 )

@test zeros(2) ∈ P2
@test length(P2) == 20
@test all( maximum( hcat(P2...), dims = 2 ) .<= 1.0 )
@test all( minimum( hcat(P2...), dims = 2 ) .>= 0.0 )

# Designs in 3D using seeds

seed_set = [ rand(3) for i = 1 : 10 ];
seed_base = deepcopy(seed_set);

# Less points requested than provided as seeds
P3_1 = monte_carlo_th( 5, 3; seeds = seed_set )
@test P3_1 == seed_set
@test seed_set == seed_base

# test whether seeds are contained in returned set
P3_2 = monte_carlo_th( 15, 3; seeds = seed_set )
#@test seed_set ⊆ P3_2
@test all( P3_2[1:10] .≈ seed_set )
@test seed_set == seed_base

# Weed out seeds that are violating box constraints
lb = .45 .* ones(3);
ub = .55 .* ones(3);

P3_3 = monte_carlo_th( 20, lb, ub; seeds = seed_set, clean_seeds = true )

@test seed_set == seed_base
@test all( maximum( hcat(P3_3...), dims = 2 ) .<= ub )
@test all( minimum( hcat(P3_3...), dims = 2 ) .>= lb )

# Keep seeds that are violating box constraints but sample only in box
P3_4 = monte_carlo_th( 20, lb, ub; seeds = seed_set, clean_seeds = false )

#@test seed_set ⊆ P3_4
@test all( P3_4[1:10] .≈ seed_set )
P3_4 = P3_4[11 : end ]
@test all( maximum( hcat(P3_3...), dims = 2 ) .<= ub )
@test all( minimum( hcat(P3_3...), dims = 2 ) .>= lb )

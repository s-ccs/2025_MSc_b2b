using Test
using Statistics
using MScB2B

@testset "Simulation sanity checks" begin

    cfg = SimulationConfig()

    epoched = simulate_cases(cfg; continuous = false, seed = 12,)   

    continuous = simulate_cases(cfg; continuous = true, seed = 12,)

    # basic dimensions
    @test size(epoched.clean[1], 1) == cfg.n_channels
    @test size(epoched.clean[1], 3) == cfg.n_trials
    @test size(continuous.clean[1], 1) == cfg.n_channels

    # clean case: predictors should be uncorrelated
    clean_evts = epoched.clean[2]
    @test abs(cor(clean_evts.condition_num, clean_evts.continuous,)) < 0.15

    # confounded case: predictors should be correlated
    confounded_evts = epoched.confounded[2]
    @test isapprox(cor(confounded_evts.condition_num, confounded_evts.continuous,), cfg.rho; atol = 0.05,)
end
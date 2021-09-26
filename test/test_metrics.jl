"""
Some tests check that KGE or KGE' is ≈ -0.41
when Q_sim is equal to μ_obs. See Knoben et al., (2019).

In other words, values > -0.41 indicates that the model
performs "better" than simply taking the mean of observations.


References
----------
1. Knoben, W.J.M., Freer, J.E., Woods, R.A., 2019. 
    Technical note: Inherent benchmark or not? 
    Comparing Nash-Sutcliffe and Kling-Gupta efficiency scores (preprint). 
    Catchment hydrology/Modelling approaches. 
    https://doi.org/10.5194/hess-2019-327
"""

using Statistics
using Test
using Streamfall


test_seq = [1,10,100,5,5,4,10,99]
b_seq = repeat([mean(test_seq)], length(test_seq))
c_seq = repeat([0.0], length(test_seq))


@testset "Kling-Gupta" begin
    @test Streamfall.KGE(test_seq, test_seq) == 1.0
    @test Streamfall.KGE(test_seq, b_seq) ≈ -0.4142135623730
    @test Streamfall.KGE(test_seq, c_seq) < 0.0
end


@testset "Modified Kling-Gupta" begin
    @test Streamfall.mKGE(test_seq, test_seq) == 1.0
    @test Streamfall.mKGE(test_seq, b_seq) ≈ -0.4142135623730
    @test Streamfall.mKGE(test_seq, c_seq) < 0.0

    
end


@testset "non-parametric Kling-Gupta" begin
    @test Streamfall.npKGE(test_seq, test_seq) == 1.0
    @test Streamfall.npKGE([10.0, 10.0, 10.0], [10.0, 10.0, 10.0]) == 1.0
    @test Streamfall.npKGE([0.0, 0.0, 0.0], [10.0, 0.0, 0.0]) == 1.0
    @test Streamfall.npKGE(test_seq, c_seq) < 0.0

    @test !isnan(Streamfall.npKGE(repeat([0], 100), repeat([0.01], 100)))
end
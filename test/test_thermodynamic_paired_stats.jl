# Run: julia --project=Analysis Analysis/test/test_thermodynamic_paired_stats.jl
using Test
include(joinpath(@__DIR__, "..", "src", "figures", "cross_so4_thermodynamics.jl"))

@testset "Exact sign test and familywise adjustment" begin
  @test direction_counts([1.0,2,3,4]).p_sign ≈ 0.125
  @test direction_counts([-1.0,1]).p_sign == 1
  @test direction_counts([0.0,0]).p_sign == 1
  @test direction_counts([0.0,1,2,3,4]).p_sign ≈ 0.125
  @test isnan(direction_counts([NaN]).p_sign)
  @test direction_counts([-0.01,0,0.01,0.2]; tolerance=0.01).n_tie == 3
  @test holm_adjust([0.01,0.04,0.03]) ≈ [0.03,0.06,0.06]
  @test holm_adjust([NaN,0.01])[2] ≈ 0.02
  @test_throws ErrorException direction_counts([1.0]; tolerance=-1)
end

@testset "Pair by seed, not by row order; preserve undefined pairs" begin
  data = DataFrame(seed=[2,1,2,1], SO4_mgL=[1000.,10,10,1000],
    category=fill("beta_oxidation",4), regime=fill("SYNTROPHY_ACT",4),
    inhibition=[0.8,0.2,0.3,0.1])
  pairs = paired_inhibition_table(data)
  beta = pairs[pairs.category .== "beta_oxidation", :]
  @test beta.seed == [1,2]
  @test beta.delta_inhibition ≈ [-0.1,0.5]
  @test count(pairs.valid_pair) == 2
  @test_throws ErrorException paired_inhibition_table(vcat(data,data[1:1,:]))
  data.inhibition[1] = NaN
  @test count(paired_inhibition_table(data).valid_pair) == 1
  draws = bootstrap_medians([0.3,0.3,0.3],MersenneTwister(1); replicates=100)
  @test all(==(0.3),draws)
  @test bootstrap_medians([-1.,0,1],MersenneTwister(7);replicates=100) ==
    bootstrap_medians([-1.,0,1],MersenneTwister(7);replicates=100)
end

@testset "Repeated-measures rank statistic with ties" begin
  data = DataFrame(seed=repeat(1:3;inner=3), SO4_mgL=repeat([10.,100,1000],3),
    category=fill("beta_oxidation",9), regime=fill("SYNTROPHY_ACT",9),
    inhibition=repeat([0.1,0.2,0.3],3))
  result = friedman_inhibition_summary(data)
  @test result.statistic[1] ≈ 6
  @test result.kendall_w[1] ≈ 1
  data.inhibition .= repeat([0.1,0.1,0.3],3)
  @test friedman_inhibition_summary(data).statistic[1] ≈ 6
  data.inhibition .= 0.1
  @test friedman_inhibition_summary(data).p_friedman_asymptotic[1] == 1
  data.inhibition[1] = NaN
  @test friedman_inhibition_summary(data).n_complete[1] == 2
end

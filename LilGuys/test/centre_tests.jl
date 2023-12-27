using LilGuys


@testset "centroid" begin
    x = [1. 5;
         0.5 1.5;
         -1 1;]

    expected = [3., 1, 0]
    cen = LilGuys.centroid(x)

    @test cen ≈ expected
end


@testset "centroid weights" begin
    x = [1. 4 10;
         0 3 0;
         -1 1 -11.3;]
    w = [1., 2, 0]

    expected = [3., 2., 1/3]

    cen = LilGuys.centroid(x, w)
    @test cen ≈ expected
end


@testset "centroid err" begin


end


@testset "centroid weights err" begin


end


@testset "centroid snap" begin
    N = 10
    pos = transpose(hcat(ones(N), zeros(N), -ones(N)))
    vel = transpose(hcat(zeros(N), 2*ones(N), zeros(N)))
    m = LilGuys.ConstVector(1., N)
    snap = LilGuys.Snapshot(pos=pos, vel=vel, m=m)

    cen = LilGuys.centroid(snap)

    @test cen.pos ≈ [1., 0., -1.]
    @test cen.vel ≈ [0., 2., 0.]
end


@testset "normal cdf" begin

    N = 1000
    μ = randn(N)
    σ = 0.5 .+ abs.(randn(N))
    x = randn(N)

    actual = LilGuys.normal_cdf.(μ, μ, σ)
    expected = fill(0.5, N)
    @test actual ≈ expected

    actual = @. LilGuys.normal_cdf(μ + σ, μ, σ) - LilGuys.normal_cdf(μ - σ, μ, σ)
    expected = fill(0.682689492137, N)
    @test actual ≈ expected

    actual = @. LilGuys.normal_cdf(μ + 100σ, μ, σ)
    expected = ones(N)
    @test actual ≈ expected

    actual = @. LilGuys.normal_cdf(μ - 100σ, μ, σ)
    expected = zeros(N)
    @test actual ≈ expected

    actual = LilGuys.normal_cdf.(μ .+ x, μ, σ)
    expected = @. 1 - LilGuys.normal_cdf(μ - x, μ, σ)
    @test actual ≈ expected
end




@testset "phase volumes" begin


end



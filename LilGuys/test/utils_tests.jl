@testset "mean" begin
    @test lguys.mean([1,2,3]) == 2
    @test lguys.mean([1, 0.5, 0, -0.5, -1]) == 0.
    @test lguys.mean([-0.3]) == -0.3
end

@testset "variance" begin
    @test lguys.var([1,2,3]) == 1 skip=true
    @test lguys.var([1, 0.5, 0, -0.5, -1]) == 0.5 * 5/4 skip=true# sample var
    @test lguys.var(fill(1, 10)) == 0 skip=true
end


@testset "percentile" begin
    x = lguys.randu(0, 1, 1000)

    qs = 100 * [0.1, 0.25, 0.5, 0.75, 0.9]
    p = lguys.percentile(x, qs)
    @test p ≈ qs/100 atol=0.1
end

@testset "randu" begin
    s = (1243, 2)
    low = -0.5
    high = 0.7
    xs = lguys.randu(low, high, s...)

    @test size(xs) == s
    μ = lguys.mean(xs)
    σ = lguys.std(xs)

    @test μ ≈ (low + high) / 2 atol=0.03
    @test σ ≈ (high - low) / √12 atol=0.03

end


@testset "rand unit vector" begin
    N = 1000
    xs = lguys.rand_unit(N)

    @test size(xs) == (3, N)
    rs = lguys.calc_r(xs)
    @test rs ≈ fill(1, N)
    μ = lguys.mean(xs)
    σ = lguys.std(xs)

    @test μ ≈ 0 atol=0.03
    @test σ > 0.3

end


@testset "gradient" begin
    x = LinRange(-1, 1, 100)
    y = x .^ 2
    g = lguys.gradient(y, x)
    @test g ≈ 2x atol=0.05

end


@testset "lerp" begin 
end
    

@testset "normal cdf" begin
    N = 1000
    μ = randn(N)
    σ = 0.5 .+ abs.(randn(N))
    x = randn(N)

    actual = lguys.normal_cdf.(μ, μ, σ)
    expected = fill(0.5, N)
    @test actual ≈ expected

    actual = @. lguys.normal_cdf(μ + σ, μ, σ) - lguys.normal_cdf(μ - σ, μ, σ)
    expected = fill(0.682689492137, N)
    @test actual ≈ expected

    actual = @. lguys.normal_cdf(μ + 100σ, μ, σ)
    expected = ones(N)
    @test actual ≈ expected

    actual = @. lguys.normal_cdf(μ - 100σ, μ, σ)
    expected = zeros(N)
    @test actual ≈ expected

    actual = lguys.normal_cdf.(μ .+ x, μ, σ)
    expected = @. 1 - lguys.normal_cdf(μ - x, μ, σ)
    @test actual ≈ expected
end


@testset "gaussian" begin
    N = 1000
    μ = randn(N)
    σ = 0.5 .+ abs.(randn(N))
    x = randn(N)
    
    actual = lguys.gaussian.(μ, μ, σ)
    expected = 1 ./ (σ * √(2π))
    @test actual ≈ expected

    actual = lguys.gaussian.(μ .- x, μ, σ) 
    expected = lguys.gaussian.(μ .+ x, μ, σ)
    @test actual ≈ expected

    actual = lguys.gaussian.(μ .+ 100σ, μ, σ)
    expected = zeros(N)
    @test actual ≈ expected

    actual = lguys.gaussian.(μ .- 100σ, μ, σ)
    expected = zeros(N)
    @test actual ≈ expected


    actual = lguys.gaussian.(μ .+ σ, μ, σ)
    expected = exp(-0.5) ./ (σ * √(2π))
    @test actual ≈ expected
end


@testset "logistic" begin
    @test lguys.logistic(0) ≈ 0.5
    @test lguys.logistic(20) ≈ 1 rtol=1e-6
    @test lguys.logistic(-20) ≈ 0 atol=1e-6

    x = LinRange(-10, 10, 100)
    y = lguys.logistic.(x)
    @test all(y .≥ 0)

end


@testset "make_equal_number_bins" begin
    x = randn(1000)
    edges = lguys.make_equal_number_bins(x, 10)
    @test length(edges) == 11
    @test issorted(edges)
    @test edges[1] == minimum(x)
    @test edges[end] == maximum(x)
end


@testset "calc_histogram" begin
    x = randn(1000)
    edges, counts = lguys.calc_histogram(x, 10)
    @test sum(counts) == length(x) - sum(x .>= maximum(x))
    @test all(counts .≥ 0)
end

    
# centroid tests
#
@testset "centroid" begin
    x = [1. 5;
         0.5 1.5;
         -1 1;]

    expected = [3., 1, 0]
    exp_err = sqrt(2) * sqrt(1/3*(2^2 + 0.5^2 + 1)) / √2 # factor for duplication and mean
    cen = lguys.centroid(x)
    err = lguys.centroid_err(x)

    @test cen ≈ expected
    @test err ≈ exp_err

    x = [1 1 1;
         2 1 0;
         -2 0 2]

    expected = [1, 1, 0]
    exp_err = sqrt(0 + 2*1^2 + 2*2^2) / sqrt(3*3) / sqrt(2)
    cen = lguys.centroid(x)
    err = lguys.centroid_err(x)

    @test cen ≈ expected
    @test err ≈ exp_err

end


@testset "centroid weights" begin
    x = [ 1. 4  10;
          0  3   0;
         -1  1 -11.3;]
    w = [1., 2,  0]

    expected = [3., 
                2., 
                1/3]

    exp_err = sqrt(1*2^2 + 2*1^2 + 
                   1*2^2 + 2*1^2 + 
                   1*(4/3)^2  + 2*(2/3)^2) / sqrt(3*3) / sqrt(2) 
    # factors for weights, variance, and sample error
    cen = lguys.centroid(x, w)
    err = lguys.centroid_err(x, w)
    @test cen ≈ expected
    @test err ≈ exp_err

    # does this reduce when weights are all equal?
    w = fill(1.23, 3)
    cen = lguys.centroid(x, w)
    err = lguys.centroid_err(x, w)
    cen2 = lguys.centroid(x)
    err2 = lguys.centroid_err(x)
    @test cen ≈ cen2
    @test err ≈ err2
end


@testset "centroid stat" begin
    N = 100000
    x = randn(3, N)
    cen = lguys.centroid(x)
    err = lguys.centroid_err(x)

    σ = 1/sqrt(N)
    @test cen ≈ zeros(3) atol=5σ
    @test err ≈ 1/sqrt(N) atol=5σ
end

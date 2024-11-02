@testset "const vector" begin
    @testset "initialization" begin
        a = 234
        N=10
        v = lguys.ConstVector(a, N)

        @test length(v) == N
        @test size(v) == (N,)

        @test eltype(v) == Float64
        @test all(v .== a)
        @test collect(v) == [a for _ in 1:N]

        #@test_throws InexactError lguys.ConstVector(a, -32)

        @test false broken = true #println
    end

    @testset "left-multiplication" begin
        a = π
        N = 4
        v = lguys.ConstVector(a, N)

        for m in [2, 0.23, exp(1)]
            @test m*a == (m*v).value
            @test m*a == (v*m).value
        end
        
    end
end



@testset "struct_to_dict" begin
    struct Foo
        a::Int
        b::Float64
        c::String
    end

    f = Foo(1, 2.0, "3")
    d = lguys.struct_to_dict(f)
    @test d == Dict(:a => 1, :b => 2.0, :c => "3")
end


@testset "dict_to_tuple" begin
    d = Dict("a" => 1, "b" => 2.0, "q" => "3")
    t = lguys.dict_to_tuple(d)

    # may be reordered, so test for set equality
    @test Set(keys(t))  == Set([:a, :b, :q])
    @test t.a == 1
    @test t.b == 2.0
    @test t.q == "3"

end

@testset "mean" begin
    # mean is imported from statsbase
    @test lguys.mean([1,2,3]) == 2
    @test lguys.mean([1, 0.5, 0, -0.5, -1]) == 0.
    @test lguys.mean([-0.3]) == -0.3
end


@testset "quantile" begin
    # percentile is imported from statsbase
    x = lguys.randu(0, 1, 1000)

    qs = [0.1, 0.25, 0.5, 0.75, 0.9]
    p = lguys.quantile(x, qs)
    @test p ≈ qs atol=0.1
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
    @test maximum(xs) ≈ high atol=0.01
    @test minimum(xs) ≈ low atol=0.01

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


@testset "gradient: simple cases" begin
    x = LinRange(-1, 1, 100)
    y = x .^ 2
    g = lguys.gradient(y, x)
    @test g ≈ 2x atol=0.05

    dx = x[2] - x[1]

    g = lguys.gradient(y)
    @test g ≈ 2x * dx atol=0.05
end





@testset "lerp" begin 
    x = [0., 0.5, 1.]
    y = [0., 1., 0.5]

    f = lguys.lerp(x, y)
    @test f.([-25, 0, 0.1, 0.3, 0.6, 1., 1.1]) ≈ [0., 0., 0.2, 0.6, 0.9, 0.5, 0.5]
    @test f(Inf) === NaN
    @test f(-Inf) === NaN
    @test f(NaN) === NaN
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



@testset "sample_Σ" begin
    @testset "uniform" begin
        low = 0.1
        high = 1

        f(x) = low < x < high ? 1/x : 0

        N = 1000

        r = lguys.sample_Σ(f, N)
        @test minimum(r) ≈ low atol=0.03
        @test maximum(r) ≈ high atol=0.03
        @test lguys.mean(r) ≈ (low + high) / 2 atol=0.02
        @test lguys.std(r) ≈ (high - low) / sqrt(12) atol=0.05
    end


    @testset "constant" begin
        f(x) = 1
        N = 1000
        log_r = log10.(LinRange(0, 1, N))
        r = lguys.sample_Σ(f, N, log_r=log_r)
        @test minimum(r) ≈ 0 atol=0.03
        @test maximum(r) ≈ 1 atol=0.03
        @test lguys.mean(r) ≈ 2/3 atol=0.03
        @test lguys.std(r) ≈ 1/4 atol=0.03
    end

end


@testset "effective_sample_size" begin
    N = 33
    w = fill(0.22, N)
    @test lguys.effective_sample_size(w) ≈ N

    w = [π; zeros(N-1)]
    @test lguys.effective_sample_size(w) ≈ 1

    w = [1,2,3]
    @test lguys.effective_sample_size(w) ≈ 36 / (1 + 2^2 + 3^2)
end

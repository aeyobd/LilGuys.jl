# these tests are written mostly for basic functionality
# and to ensure the interface is as documented.
# Correctness is assumed since these are functions 
# from other packages (

@testset "stats" begin
    x = [1, -0.1, 0.3]
    w = [0.6, 0, 0.3]

    @test lguys.mean(x) ≈ 0.4
    @test lguys.mean(x, w) ≈ 0.69 / 0.9
    
    expected = sqrt(((1-0.4)^2 + (-0.1-0.4)^2 + (0.3-0.4)^2)/2)
    @test lguys.std(x) ≈ expected
    @test lguys.variance(x) ≈ expected^2

    s = sqrt( ((1-0.7666666666666666)^2 * 0.6 +(0.3 - 0.7666666666666666)^2 * 0.3) / 0.9)
    @test lguys.std(x, w) ≈ s
    @test lguys.variance(x, w) ≈s^2

    @test lguys.quantile(x, 0.5) ≈ 0.3
    @test lguys.quantile(x, 0.) ≈ -0.1
    @test lguys.quantile(x, w, 1.) ≈ 1
end

@testset "midpoints" begin
    x = [1,2,3]
    @test LilGuys.midpoints(x) ≈ [1.5, 2.5]

    x = [-2, 0., 2., 2.5]
    @test LilGuys.midpoints(x) ≈ [ -1, 1, 2.25]
end


@testset "expint" begin
    @test LilGuys.expinti(0) === -Inf
    @test LilGuys.expinti(1) ≈ 1.89511781635594 # sage
    @test LilGuys.expinti(3) ≈ 9.93383257062542
end


@testset "erf" begin
    @test lguys.erf(1) ≈ 0.842700792949715 # sage
    @test lguys.erf(-2.5) ≈ -0.999593047982555
    @test lguys.erf(0) ≈ 0.0
end


@testset "integrate" begin
    let
        f(x) = x^2
        @test lguys.integrate(f, 0, 1) ≈ 1/3
    end

    let
        f(x) = 1/x
        @test lguys.integrate(f, 1e-30, 1) > 0
    end
end


@testset "curve_fit" begin
    let 
        f(x, p0) = @. p0[1] + p0[2] * x + p0[3] * x^2
        popt, covt = lguys.curve_fit(f, [1., 2., 3., 4.], [1., 2., 4., 7.], [0., 0., 0.])
        @test popt ≈ [1., -1/2, 1/2]
    end

    let 
        f(x, p0) = @. p0[1] + p0[2] * x + p0[3] * x^2
        popt, covt = lguys.curve_fit(f, [1., 2., 3., 4.], [1., 2., 4., -12.], [1, 0.5, 0.2, 0], [0., 0., 0.])
        @test popt ≈ [1., -1/2, 1/2]
    end
end


@testset "filter_nans" begin
    x = [1, 2, NaN, 3, Inf]
    y, w = lguys.Interface.filter_nans(x)
    @test y ≈ [1., 2, 3]

    x = [1, 2, NaN, 3, Inf]
    w = [0.5, 0.5, 0.5, NaN, 0.5]
    y, w = lguys.Interface.filter_nans(x, w)
    @test y ≈ [1., 2]
    @test w ≈ [0.5, 0.5]
end



@testset "histogram" begin
    x = [1,1.5, 2, 3]
    bins = [0, 2.5, 5]
    bins, values, err = lguys.histogram(x, bins)
    @test bins ≈ [0, 2.5, 5]
    @test values ≈ [3, 1]
    @test err ≈ [√3, 1]
end



@testset "bins_default" begin
    x = [1,1.5, 2, 3, 4]
    w = 2 * (1.5) / cbrt(5)

    bins =  lguys.Interface.bins_default(x, nothing)
    l = ceil((4-1 + w)/w)
    @test length(bins) == l + 1
    @test bins[1] < 1
    @test bins[end] > 4
    @test bins[2] - bins[1] ≈ w
end

@testset "bins_equal_width" begin
    x = [3,1,7]
    w = 2.
    bins = lguys.Interface.bins_equal_width(x, nothing; bin_width=w)
    @test bins ≈ [0, 2, 4, 6, 8]

end


@testset "bins_equal_number" begin
    x = [3, 1, 7]
    bins = lguys.Interface.bins_equal_number(x, nothing; num_per_bin=1.5)
    @test bins ≈ [1, 3, 7]
end


@testset "bins_both" begin
    x = randn(100)

    bins = lguys.Interface.bins_both(x, nothing; num_per_bin=9, bin_width=0.5)
    @test length(bins) > 1
    
    _, values, _ = lguys.histogram(x, bins)

    @test sum(values) == 100
    @test all(values .>= 9)
    @test all(diff(bins) .>= 0.5)

end


@testset "default_bin_width" begin
    x = [1,1.5, 2, 3, 4]
    @test lguys.default_bin_width(x) ≈ 2 * (1.5) / cbrt(5)
end


@testset "default_n_per_bin" begin
    x = randn(100)
    @test lguys.default_n_per_bin(x) ≈ (2*100^(2/5))
end





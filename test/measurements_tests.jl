@testset "initialization" begin
    meas = Measurement(1, 0.2)
    @test meas isa Measurement{<:Float64}
    @test meas.middle == 1
    @test meas.lower == 0.2
    @test meas.upper == 0.2
    @test meas.kind == ""
    @test lguys.sym_error(meas) == 0.2

    meas = Measurement{Float64}(1, 0, 2, "test")
    @test meas.middle == 1
    @test meas.lower == 0
    @test meas.upper == 2
    @test meas.kind == "test"
    @test lguys.sym_error(meas) == 2


    meas = Measurement(0.5)
    @test meas.middle == 0.5
    @test error_interval(meas) == (0.0, 0.0)

end


@testset "addition" begin
    @testset "special cases" begin
        a = Measurement(10, 1.1, 0.9)
        b = Measurement(1.8, 2.2, "normal")
        @test_logs (:warn, "measurement kinds do not match") a+b

        a = Measurement(10, 1.1, 0.9, "normal")
        @test a + b ≈ Measurement(11.8, 3.3, 3.1)
        @test 0.5 + a ≈ Measurement(10.5, 1.1, 0.9)
        @test (a+b).kind == "normal"
        @test (a+1).kind == "normal"

        b = Measurement(1.8, 2.2, "normal")
        @test a - 2.1 ≈ Measurement(7.9, 1.1, 0.9)
        @test a - b ≈ Measurement(8.2, 3.3, 3.1)
    end

end


@testset "isapprox" begin
    a = Measurement(10., 1.1)
    b = Measurement(10., 1.1 * (1+eps()))
    @test a ≈ b

    a = Measurement(10., 1.1)
    b = Measurement(10., 1.1 * (1+cbrt(eps())))
    @test !(a ≈ b)

    a = Measurement(10., 1.1)
    b = Measurement(9.9984, 1.1)
    @test !(a ≈ b)
    @test a ≈ b atol=2e-3
    @test a ≈ b rtol=2e-2
    @test b ≈ 9.9984
    @test !(b ≈ 9.8)


    a = Measurement(3.2e9, NaN, 1e9)
    b = Measurement(3.1e9, NaN, 1.2e9)
    @test a ≈ b nans=true rtol=0.2

    a = Measurement(3.2e9, -Inf, NaN)
    b = Measurement(3.1e9, -Inf, NaN)
    @test a ≈ b nans=true rtol=0.2
    @test !isapprox(a, b, rtol=0.2)
    @test a ≈ 3.2e9

    a = Measurement(NaN, 0.1, NaN)
    b = a
    @test a ≈ b nans=true
    @test !isapprox(a, b)
    @test a ≈ NaN nans=true


    a = Measurement(104.2, 0.2, 0.4)
    b = Measurement(104.21, 0.19, 0.38)
    @test a ≈ b atol = 0.03
    @test a ≈ b rtol = 0.03 broken=true # would like to implement rtol relative to middle :/
    @test !(a≈b)
end



@testset "log10" begin
    @testset "special cases" begin
        a = Measurement(20., 1.5)
        b = log10(a)
        @test b ≈ Measurement(1.3010, 0.034, 0.031) atol=0.001

        a = Measurement(3)
        b = Measurement(log10(3))
        @test log10(a) ≈ b

        a = Measurement(1.2, 1.2, 1.5)
        b = Measurement(0.07918, Inf, 0.35218) 
        @test log10(a) ≈ b atol=1e-4
    end
end


@testset "logexp inverse" begin
    for i in 1:100
        a = Measurement(3*randn(), 0.2*rand(), 0.2*rand())
        b = log10(exp10(a))
        @test a ≈ b
    end
end


@testset "sqrt sq inverse" begin
    N = 100
    x0 = 10 .^ (4.5 * randn(N))
    a = Measurement.(x0, rand(N) .*x0 / 5, rand(N) .*x0/8)
    b = @. sqrt(a^2)
    @test a ≈ b
    @test all(b .≈ x0)
end


@testset "exponentiation" begin
    a = Measurement(0.3, 0.1)
    b = Measurement(1.995262314, 0.410369, 0.516624)

    @test b ≈ 10^a atol = 1e-5

    b = Measurement(1.31638, 0.11526, 0.1263177)
    @test 2.5 ^ a ≈ b atol = 1e-5
end


@testset "isfinite" begin
    a = Measurement(NaN, 0.1, 0.3)
    @test !isfinite(a)
    a = Measurement(2.3, 0.1)
    @test isfinite(a)

    a = Measurement(0.5, -Inf, 2)
    @test !isfinite(a)

    a = Measurement(2.2, 0, -NaN)
    @test !isfinite(a)
end


@testset "comparisons" begin
    a = Measurement(2.3, 0.1)
    b = Measurement(3.0, 2, π)
    @test a < b
    @test -a > -b
    @test a < 3
    @test 2.99 <= b
    @test !(NaN < a)
end

@testset "negative" begin
    a = Measurement(-105, 0.5, 0.8)
    @test -a ≈ Measurement(105, 0.5, 0.8)
    @test -a ≈ 0 - a
end

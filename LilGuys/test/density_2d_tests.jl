function calc_χ2(x, x_exp, xerr)
    filt = @. isfinite(x) & isfinite(x_exp) & isfinite(xerr)
	return sum(
               @. ((x-x_exp)^2 / xerr^2)[filt]
        ) / sum(filt)
end



@testset "calc_Σ" begin
    log_r = [-Inf, 0]
    mass = [1]
    Σ = LilGuys.calc_Σ_from_hist(log_r, mass)
    @test Σ ≈ [1/π]

    log_r = [-1, 0, 1, 1.5]
    mass = [1, 2, 0]
    Σ = LilGuys.calc_Σ_from_hist(log_r, mass)
    @test Σ ≈ [1 / (0.99π), 2 / (99π), 0]
end


@testset "calc_Σ_from_1D_density" begin
end

@testset "calc_Σ_mean" begin
    log_r = [-Inf, 0]
    mass = [1]
    Σ = LilGuys.calc_Σ_mean_from_hist(log_r, mass)
    @test Σ ≈ [1/π]

    log_r = [-1, 0, 1, 1.5]
    mass = [1, 2, 0]
    mass = cumsum(mass)
    Σ = LilGuys.calc_Σ_mean_from_hist(log_r, mass)
    @test Σ ≈ [1 / (π), 3 / (100π), 3 / (1000π)]
end


@testset "integration with exp profile" begin
    Σ(r) = exp(-r)/2π
    N = 10_000
    r = LilGuys.sample_Σ(Σ, N, log_r=LinRange(-5, 5, 1000))

    mass = 0.5 .+ 0.5rand(N)
    M = sum(mass)

    obs = LilGuys.calc_properties(r, normalization=:none, weights=mass)

    @test sum(obs.mass_in_annulus) ≈ sum(mass)
    @test obs.M_in[end] ≈ sum(mass)
    @test issorted(obs.log_r)
    @test sum(obs.counts) ≈ N

    r = 10 .^ obs.log_r
    sigma_exp = M * Σ.(r)
    χ2 = calc_χ2(obs.Sigma, sigma_exp, obs.Sigma_err)
    println(obs.Sigma)
    println(sigma_exp)
    println(obs.Sigma_err)
    @test χ2  ≈ 1 rtol = 0.5

    Gamma_exp = -r
    χ2 = calc_χ2(obs.Gamma, Gamma_exp, obs.Gamma_err)
    @test χ2  ≈ 1 rtol = 0.5


    r = 10 .^ obs.log_r_bins[2:end]
    M_exp = @. M * (1 - exp(-r) - r*exp(-r))
    χ2 = calc_χ2(obs.M_in, M_exp, obs.M_in_err)
    @test χ2  ≈ 1 rtol = 0.7


    @testset "read/write" begin
        dir = mktempdir()

        filename = joinpath(dir, "test.toml")
        open(filename, "w") do f
            LilGuys.print(f, obs)
        end

        obs2 = LilGuys.ObsProfile(filename)

        for name in fieldnames(typeof(obs))
            v = getproperty(obs, name)
            if v isa String
                @test getproperty(obs2, name) == v
            else    
                filt = .!isnan.(v)
                @test v[filt] ≈ getproperty(obs2, name)[filt]
                @test isnan.(v) == isnan.(getproperty(obs2, name))
            end
        end
    end
end



@testset "calc_r_ell" begin
    N = 100
    @testset "circular" begin
        x = randn(N)
        y = randn(N)
        r = sqrt.(x.^2 + y.^2)
        r_ell = LilGuys.calc_r_ell(x, y, 0, 23.425)

        @test r_ell ≈ r

        a = rand(N)
        b = a
        PA = 360rand(N)
        r_ell = LilGuys.calc_r_ell.(x, y, a, b, PA)
        @test r_ell ≈ r ./ a
    end

    @testset "PA = 0, 90" begin
        x = randn(N)
        y = randn(N)
        a = rand(N)
        b = rand(N)

        # PA of zero points north, so major axis is y
        r = @. sqrt(x^2/b^2 + y^2/a^2)
        r_ell = LilGuys.calc_r_ell.(x, y, a, b, 0.)
        @test r_ell ≈ r

        # PA of 90 points east, so major axis is x
        r = @. sqrt(x^2/a^2 + y^2/b^2)
        r_ell = LilGuys.calc_r_ell.(x, y, a, b, 90.)
        @test r_ell ≈ r
    end

    @testset "zero" begin
        x = 0
        y = 0
        ell = rand(N)
        PA = 360rand(N)

        r = zeros(N)
        r_ell = LilGuys.calc_r_ell.(x, y, ell, PA)
        @test r_ell ≈ r

        a = rand(N)
        b = rand(N)
        r_ell = LilGuys.calc_r_ell.(x, y, a, b, PA)
        @test r_ell ≈ r
    end

    @testset "periodicic PA" begin
        x = randn(N)
        y = randn(N)
        a = rand(N)
        b = rand(N)
        PA = 360rand(N)

        r1 = LilGuys.calc_r_ell.(x, y, a, b, PA)
        r2 = LilGuys.calc_r_ell.(x, y, a, b, PA .+ 360)
        r2 = LilGuys.calc_r_ell.(x, y, a, b, PA .+ 180)
        r3 = LilGuys.calc_r_ell.(x, y, a, b, PA .- 360)
        r4 = LilGuys.calc_r_ell.(x, y, a, b, PA .+ 720)
        @test r1 ≈ r2
        @test r1 ≈ r3
        @test r1 ≈ r4
    end

    @testset "exceptions" begin
        @test_throws DomainError LilGuys.calc_r_ell(1., 2, 0, 0, 90)
        @test_throws DomainError LilGuys.calc_r_ell(1., 2, 0, 1, 90)
        @test_throws DomainError LilGuys.calc_r_ell(1., 2, 1, 0, 90)
        @test_throws DomainError LilGuys.calc_r_ell(1., 2, 1, -1, 90)
    end
end


@testset "calc_r_ell_sky" begin
    @test false broken=true
end

@testset "shear_points_to_ellipse" begin
    @test false broken=true
end


@testset "spherical_mean" begin
    @test false broken=true
end


@testset "calc_centre2D" begin
    @test false broken=true
end

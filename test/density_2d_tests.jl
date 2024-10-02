

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

    obs = LilGuys.StellarProfile(r, normalization=:none, weights=mass)

    @test sum(obs.mass_in_annulus) ≈ sum(mass)
    @test obs.M_in[end] ≈ sum(mass)
    @test issorted(obs.log_r)
    @test sum(obs.counts) ≈ N

    r = 10 .^ obs.log_r
    sigma_exp = M * Σ.(r)
    @test_χ2 obs.Sigma obs.Sigma_err sigma_exp

    Gamma_exp = -r
    @test_χ2 obs.Gamma obs.Gamma_err Gamma_exp

    r = 10 .^ obs.log_r_bins[2:end]
    M_exp = @. M * (1 - exp(-r) - r*exp(-r))
    @test_χ2 obs.M_in obs.M_in_err M_exp


    @testset "read/write" begin
        dir = mktempdir()

        filename = joinpath(dir, "test.toml")
        open(filename, "w") do f
            LilGuys.print(f, obs)
        end

        obs2 = LilGuys.StellarProfile(filename)

        for name in fieldnames(typeof(obs))
            v = getproperty(obs, name)
            if v isa String
                @test getproperty(obs2, name) == v
            elseif v isa Real
                if isnan(v)
                    @test isnan(getproperty(obs2, name))
                else
                    @test v ≈ getproperty(obs2, name)
                end
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


@testset "to_orbit_coords" begin
    @testset "identity" begin
        N = 1000
        ra = -180 .+ 360rand(N)
        dec = 180rand(N) .- 90

        ra2, dec2 = lguys.to_orbit_coords(ra, dec, 0., 0., 90.)

        @test ra2 ≈ ra
        @test dec2 ≈ dec

        radec = lguys.to_orbit_coords.(ra, dec, ra, dec, 0.)
        ra3 = first.(radec)
        dec3 = last.(radec)
        @test ra3 ≈ zeros(N) atol = 1e-12
        @test dec3 ≈ zeros(N) atol = 1e-12
    end

    @testset "direction" begin
        N = 1000
        ra = -180 .+ 360rand(N)
        dec = 160rand(N) .- 80

        # these have to be very small for this to work
        δ = 0.01 * rand(N)
        α = 0.01 * rand(N)

        ra1 = ra .+ α
        dec1 = dec .+ δ
        α_exp = lguys.angular_distance.(ra, dec, ra1, dec)
        δ_exp = lguys.angular_distance.(ra, dec, ra, dec1)

        # if the PA is 0, then positive in the new axis should be positive in dec
        radec = lguys.to_orbit_coords.(ra, dec1, ra, dec, 0.)
        ra2 = first.(radec)
        dec2 = last.(radec)
        @test ra2 ≈ δ_exp atol=1e-6
        @test dec2 ≈ zeros(N) atol=1e-5


        radec = lguys.to_orbit_coords.(ra1, dec, ra, dec, 0.)
        ra2 = first.(radec)
        dec2 = last.(radec)
        @test ra2 ≈ zeros(N) atol=1e-5
        @test dec2 ≈ -α_exp atol=1e-6
        

        # if PA is 90, then positive in the new axis should be positive in ra
        radec = lguys.to_orbit_coords.(ra1, dec, ra, dec, 90.)
        ra2 = first.(radec)
        dec2 = last.(radec)
        @test ra2 ≈ α_exp atol=1e-6
        @test dec2 ≈ zeros(N) atol=1e-5

        radec = lguys.to_orbit_coords.(ra, dec1, ra, dec, 90.)
        ra2 = first.(radec)
        dec2 = last.(radec)
        @test ra2 ≈ zeros(N) atol=1e-6
        @test dec2 ≈ δ_exp atol=1e-5

        # 180 deg
        radec = lguys.to_orbit_coords.(ra, dec1, ra, dec, 180)
        ra2 = first.(radec)
        dec2 = last.(radec)
        @test ra2 ≈ -δ_exp atol=1e-6
        @test dec2 ≈ zeros(N) atol=1e-5

        radec = lguys.to_orbit_coords.(ra1, dec, ra, dec, 180)
        ra2 = first.(radec)
        dec2 = last.(radec)
        @test ra2 ≈ zeros(N) atol=1e-5
        @test dec2 ≈ α_exp atol=1e-6
    end
end



@testset "calc_Σ" begin
    log_R = [-Inf, 0]
    mass = [1]
    Σ, e = LilGuys.calc_Σ_from_hist(log_R, mass, sqrt.(mass))
    @test Σ ≈ [1/π]

    log_R = [-1, 0, 1, 1.5]
    mass = [1, 2, 0]
    Σ, e = LilGuys.calc_Σ_from_hist(log_R, mass, sqrt.(mass))
    @test Σ ≈ [1 / (0.99π), 2 / (99π), 0]
end


@testset "calc_Σ_from_1D_density" begin
end

@testset "calc_Σ_mean" begin
    log_R = [-1., 0.]
    mass = cumsum([0, 1])

    Σ = LilGuys.calc_Σ_mean_from_hist(log_R, mass)
    @test Σ ≈ [0, 1/π]

    log_R = [0, 1, 1.5]
    mass = [1, 2, 0]
    mass = cumsum(mass)
    Σ = LilGuys.calc_Σ_mean_from_hist(log_R, mass)
    @test Σ ≈ [1 / (π), 3 / (100π), 3 / (1000π)]
end


@testset "stellar profile" begin
    @test_throws DomainError LilGuys.StellarDensityProfile([-1, 2, 3])
    @test_throws ArgumentError LilGuys.StellarDensityProfile([1])

    @testset "normalization" begin
        radii = [1, 3, 4, 4.5, 6, 10]
        bins = log10.([0, 3.5, 7, 11, 20])
        r_bins = 10 .^ bins
        r_m = LilGuys.midpoints(bins)
        # 2 count in first bin (area π(3.5)^2)
        # 3 counts in second bin (area πr^2 - π(3.5)^2)
        # 1 count in third bin (area π11^2 - π7^2)
        # 0 counts in fourth bin

        areas = diff((10 .^ bins) .^ 2) * π
        mass_in_annulus_unnormed = [2,3,1,0]
        counts = [2,3,1, 0]
        mass_in_annulus_err = sqrt.(counts)
        Sigma = [2,3,1,0] ./ areas
        Sigma_err = mass_in_annulus_err ./ areas

        log_Sigma = log10.(Sigma)
        log_Sigma_em = log_Sigma .- log10.(Sigma .- Sigma_err)


        # default is no normalization
        prof = LilGuys.StellarDensityProfile(radii, bins=bins)
        @test LilGuys.middle.(prof.log_Sigma) ≈ log_Sigma nans=true
        e = LilGuys.lower_bound.(prof.log_Sigma)
        @test e[1:3] ≈ log_Sigma_em[1:3] atol=1e-8

        prof = LilGuys.StellarDensityProfile(radii, bins=bins, normalization=:mass)
        Mtot = 6
        e = LilGuys.lower_bound.(prof.log_Sigma)
        @test LilGuys.middle.(prof.log_Sigma) ≈ log_Sigma .- log10(Mtot) nans=true atol=1e-8
        @test e[1:3] ≈ log_Sigma_em[1:3] atol=1e-8


        prof = LilGuys.StellarDensityProfile(radii, bins=bins, normalization=:central, bins_centre=1)
        Mtot = 2 / (π * 3.5^2)
        @test LilGuys.middle.(prof.log_Sigma) ≈ log_Sigma .- log10(Mtot) nans=true atol=1e-8
        e = LilGuys.lower_bound.(prof.log_Sigma)
        @test e[1:3] ≈ log_Sigma_em[1:3] atol=1e-8


        @test_throws Exception LilGuys.StellarDensityProfile(radii, bins=bins, normalization=:jaberwocky)
    end



    @testset "snapshot" begin
        x = [0., 0.5, 1.0, -1.0]
        y = [1.0, 0.0, -2.5, 3.0]
        z = [0.0, 1.0, 0.5, -1.0]

        positions = [x y z]'
        velocities = zeros(3, 4)
        masses = [1, 0, 1, 0]

        snap = LilGuys.Snapshot(positions, velocities, ones(4))
        snap.weights = masses

        prof = LilGuys.StellarDensityProfile(snap)

        prof = LilGuys.StellarDensityProfile(snap, x_vec = [0, 1, 0], y_vec = [0, 0, 1])

        prof = LilGuys.StellarDensityProfile(snap, R_units="arcmin")

        prof = LilGuys.StellarDensityProfile(snap, normalization=:none)
        @test false broken=true
    end
end


@testset "integration with exp profile unweighted" begin
    Σ(R) = exp(-R)/2π
    N = 10_000
    R = LilGuys.sample_Σ(Σ, N, log_R=LinRange(-5, 5, 1000))

    obs = LilGuys.StellarDensityProfile(R, normalization=:none)
    
    Σs = 10 .^ LilGuys.middle.(obs.log_Sigma)
    mass_per_annulus = sum(Σs .* diff(π * 10 .^ 2obs.log_R_bins))

    @test sum(mass_per_annulus) ≈ N
    @test issorted(obs.log_R)
    @test sum(obs.counts) ≈ N

    R = 10 .^ obs.log_R
    sigma_exp = N * Σ.(R)
    log_sigma_exp = log10.(sigma_exp)
    
    err = maximum.(LilGuys.credible_interval.(obs.log_Sigma))
    @test_χ2 LilGuys.middle.(obs.log_Sigma) err log_sigma_exp

    Gamma_exp = -R
    err = maximum.(LilGuys.credible_interval.(obs.Gamma))
    filt = isfinite.(err)
    @test sum(filt) > 10

    x2 = (LilGuys.middle.(obs.Gamma)[filt]  .- Gamma_exp[filt]) ./ err[filt]
    @info obs.Gamma[filt][argmax(x2)]
    @info err[filt][argmax(x2)]
    @info Gamma_exp[filt][argmax(x2)]

    @test_χ2 LilGuys.middle.(obs.Gamma)[filt] err[filt]  Gamma_exp[filt]


    @testset "read/write" begin
        dir = mktempdir()

        filename = joinpath(dir, "test.toml")
        open(filename, "w") do f
            LilGuys.print(f, obs)
        end

        obs2 = LilGuys.StellarDensityProfile(filename)

        for name in fieldnames(typeof(obs))
            v = getproperty(obs, name)
            if v isa String
                @test getproperty(obs2, name) == v
            elseif v isa Dict
                @test getproperty(obs2, name) == v
            elseif v isa Real
                @test v ≈  getproperty(obs2, name) nans=true
            else    
                @test all(isapprox.(v,  getproperty(obs2, name), nans=true))
            end
        end
    end
end



@testset "shear_points_to_ellipse" begin
    @testset "exceptions" begin
        @test_throws DomainError LilGuys.shear_points_to_ellipse(1, 2, 0, 0, 90)
        @test_throws DomainError LilGuys.shear_points_to_ellipse(1, 2, 0, -1, 90)
    end

    @testset "special cases" begin
        a, b = LilGuys.shear_points_to_ellipse(0, 0, √5, 0.234, 341)
        @test a ≈ 0
        @test b ≈ 0


        xs = [1, 0, -1, 0, 1]
        ys = [0, 1, 0, -1, 1]
        # major axis points +y
        # minor axis points -x
        a, b = LilGuys.shear_points_to_ellipse(xs, ys, 2, 1, 0)
        @test a ≈ ys ./ 2
        @test b ≈ -xs


        # major axis points +x
        # minor axis points +y
        a, b = LilGuys.shear_points_to_ellipse(xs, ys, 2, 1, 90)
        @test a ≈ xs ./ 2
        @test b ≈ ys


        # major axis points -y
        # minor axis points +x
        a, b = LilGuys.shear_points_to_ellipse(xs, ys, 2, 1, 180)
        @test a ≈ -ys ./ 2
        @test b ≈ xs


        # major axis points -x
        # minor axis points -y
        a, b = LilGuys.shear_points_to_ellipse(xs, ys, 2, 1, -90)
        @test a ≈ -xs ./ 2
        @test b ≈ -ys
    end

end


@testset "spherical_mean" begin
    @testset "exceptions" begin
        @test_throws DomainError LilGuys.spherical_mean([0, 24], [90, -90])
        @test_throws DimensionMismatch LilGuys.spherical_mean([0, 24], [-90])
    end

    @testset "special cases" begin
        ra, dec = LilGuys.spherical_mean([0, 0], [90, 0])
        @test ra ≈ 0
        @test dec ≈ 45

        ra, dec = LilGuys.spherical_mean([0, -40], [0, 0])
        @test ra ≈ 360 + -20
        @test dec ≈ 0

        ra, dec = LilGuys.spherical_mean([231, NaN], [23, 23])
        @test isnan(ra)
        @test isnan(dec)

        ra, dec = LilGuys.spherical_mean([231, 342], [NaN, 23])
        @test isnan(ra)
        @test isnan(dec)
    end



    @testset "random" begin 
        N = 100000

        ra0 = 23.85
        dec0 = -13.2

        ra = ra0 .+ 10randn(N)
        dec = dec0 .+ 7randn(N)

        ra_mean, dec_mean = LilGuys.spherical_mean(ra, dec)

        @test ra_mean ≈ ra0 atol=3e-1
        @test dec_mean ≈ dec0 atol=3e-1


        ra0 = 85
        dec0 = 90

        ra = 360rand(N)
        dec = dec0 .- 5rand(N)

        ra_mean, dec_mean = LilGuys.spherical_mean(ra, dec)

        @test dec_mean ≈ dec0 atol=3e-1
    end
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

@testset "integration with exp profile" begin
    Σ(r) = exp(-r)/2π
    N = 10_000
    r = LilGuys.sample_Σ(Σ, N, log_R=LinRange(-5, 5, 1000))

    mass = 0.5 .+ 0.5rand(N)
    M = sum(mass)

    obs = LilGuys.StellarDensityProfile(r, normalization=:none, weights=mass)

    mass_per_annulus = 10 .^ obs.log_Sigma .* diff(π * 10 .^ 2obs.log_R_bins)
    @test sum(LilGuys.middle.(mass_per_annulus)) ≈ sum(mass)

    @test issorted(obs.log_R)
    @test sum(obs.counts) ≈ N

    r = 10 .^ obs.log_R
    sigma_exp = M * Σ.(r)

    err = maximum.(LilGuys.credible_interval.(obs.log_Sigma))
    @test_χ2 LilGuys.middle.(obs.log_Sigma) err log10.(sigma_exp)

    Gamma_exp = -r
    err = maximum.(LilGuys.credible_interval.(obs.Gamma))
    @test_χ2 LilGuys.middle.(obs.Gamma) err  Gamma_exp



    @testset "read/write" begin
        dir = mktempdir()

        filename = joinpath(dir, "test.toml")
        open(filename, "w") do f
            LilGuys.print(f, obs)
        end

        obs2 = LilGuys.StellarDensityProfile(filename)

        for name in fieldnames(typeof(obs))
            v = getproperty(obs, name)
            if v isa String
                @test getproperty(obs2, name) == v
            elseif v isa Dict
                @test getproperty(obs2, name) == v
            elseif v isa Real
                @test v ≈ getproperty(obs2, name) nans=true
            else    
                @test all(isapprox.(v, getproperty(obs2, name), nans=true))
            end
        end
    end
end



@testset "calc_R_ell" begin
    N = 100
    @testset "circular" begin
        x = randn(N)
        y = randn(N)
        r = sqrt.(x.^2 + y.^2)
        r_ell = LilGuys.calc_R_ell(x, y, 0, 23.425)

        @test r_ell ≈ r

        a = rand(N)
        b = a
        PA = 360rand(N)
        r_ell = LilGuys.calc_R_ell.(x, y, a, b, PA)
        @test r_ell ≈ r ./ a
    end

    @testset "PA = 0, 90" begin
        x = randn(N)
        y = randn(N)
        a = rand(N)
        b = rand(N)

        # PA of zero points north, so major axis is y
        r = @. sqrt(x^2/b^2 + y^2/a^2)
        r_ell = LilGuys.calc_R_ell.(x, y, a, b, 0.)
        @test r_ell ≈ r

        # PA of 90 points east, so major axis is x
        r = @. sqrt(x^2/a^2 + y^2/b^2)
        r_ell = LilGuys.calc_R_ell.(x, y, a, b, 90.)
        @test r_ell ≈ r
    end

    @testset "zero" begin
        x = 0
        y = 0
        ell = rand(N)
        PA = 360rand(N)

        r = zeros(N)
        r_ell = LilGuys.calc_R_ell.(x, y, ell, PA)
        @test r_ell ≈ r

        a = rand(N)
        b = rand(N)
        r_ell = LilGuys.calc_R_ell.(x, y, a, b, PA)
        @test r_ell ≈ r
    end

    @testset "periodicic PA" begin
        x = randn(N)
        y = randn(N)
        a = rand(N)
        b = rand(N)
        PA = 360rand(N)

        r1 = LilGuys.calc_R_ell.(x, y, a, b, PA)
        r2 = LilGuys.calc_R_ell.(x, y, a, b, PA .+ 360)
        r2 = LilGuys.calc_R_ell.(x, y, a, b, PA .+ 180)
        r3 = LilGuys.calc_R_ell.(x, y, a, b, PA .- 360)
        r4 = LilGuys.calc_R_ell.(x, y, a, b, PA .+ 720)
        @test r1 ≈ r2
        @test r1 ≈ r3
        @test r1 ≈ r4
    end

    @testset "exceptions" begin
        @test_throws DomainError LilGuys.calc_R_ell(1., 2, 0, 0, 90)
        @test_throws DomainError LilGuys.calc_R_ell(1., 2, 0, 1, 90)
        @test_throws DomainError LilGuys.calc_R_ell(1., 2, 1, 0, 90)
        @test_throws DomainError LilGuys.calc_R_ell(1., 2, 1, -1, 90)
    end
end




@testset "calc_centre2D" begin
    ra = [23, 24]
    dec = [0, 360]
    ra_centre, dec_centre = LilGuys.calc_centre2D(ra, dec, "mean")
    @test ra_centre ≈ 23.5
    @test dec_centre ≈ 0


    weights = [1, 0]
    ra_centre, dec_centre = LilGuys.calc_centre2D(ra, dec, "mean", weights)
    @test ra_centre ≈ 23
    @test dec_centre ≈ 0

    ra_centre, dec_centre = LilGuys.calc_centre2D(ra, dec, (2, 3))
    @test ra_centre ≈ 2
    @test dec_centre ≈ 3


    @test_throws ArgumentError LilGuys.calc_centre2D(ra, dec, "jabberwocky")
end


@testset "calc_R_ell_sky" begin
    # this function integrates calc_R_ell_sky, to_tangent, and calc_centre2D
    #
    r0 = 0.01
    ra = [0., 1, 0, -1, 0] .* r0 .+ 34
    dec = [0., 0, 1, 0, -1] .* r0

    r_ell = LilGuys.calc_R_ell_sky(ra, dec, 1, 1, 90)
    @test r_ell ≈ [0., 1, 1, 1, 1] .* 60 .* r0

    r_ell = LilGuys.calc_R_ell_sky(ra, dec)
    @test r_ell ≈ [0., 1, 1, 1, 1] .* 60 .* r0

    r_ell = LilGuys.calc_R_ell_sky(ra, dec, 0, 45)
    @test r_ell ≈ [0., 1, 1, 1, 1] .* 60 .* r0


    r_ell = LilGuys.calc_R_ell_sky(ra, dec, 3, 1, 90, units="deg")
    @test r_ell ≈ [0., 1/3, 1, 1/3, 1]  .* r0

    r_ell = LilGuys.calc_R_ell_sky(ra, dec, 0.5, 0, units="deg")
    @test r_ell ≈ [0., 2, 1, 2, 1] / sqrt(2)  .* r0

    r_ell = LilGuys.calc_R_ell_sky(ra, dec, 1, 1, -24, centre=(34 - 2r0, 0), units="deg")
    @test r_ell ≈ [2., 3., √5, 1, √5]  .* r0 rtol=1e-6


    @testset "units" begin

        r_ell = LilGuys.calc_R_ell_sky(ra, dec, 1, 1, 90, units="arcmin")
        @test r_ell ≈ [0., 1, 1, 1, 1] .* 60 .* r0

        r_ell = LilGuys.calc_R_ell_sky(ra, dec, 1, 1, 23, units="arcsec")
        @test r_ell ≈ [0., 1, 1, 1, 1] .* 60 * 60 .* r0

        r_ell = LilGuys.calc_R_ell_sky(ra, dec, 1, 1, -234, units="deg")
        @test r_ell ≈ [0., 1, 1, 1, 1]  .* r0

        @test_throws Exception LilGuys.calc_R_ell_sky(ra, dec, 1, 1, 90, units="jabberwocky")
    end
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



@testset "scale" begin
    N = 100

    r = 10 .^ (0 .+ 0.4randn(N))
    m = rand(N)

    bins = LinRange(-1, 1, 5)
    prof_1 = LilGuys.StellarDensityProfile(r, weights=m, bins=bins, normalization=:none)

    M_scale = 0.232
    r_scale = 1.992
    prof_1 = LilGuys.scale(prof_1, r_scale, M_scale)

    prof_2 = LilGuys.StellarDensityProfile(r * r_scale, weights=m * M_scale, bins=bins .+ log10(r_scale), normalization=:none)

    for key in fieldnames(typeof(prof_1))
        v1 = getproperty(prof_1, key)
        v2 = getproperty(prof_2, key)
        if key ∈ [:log_R_scale, :log_m_scale]
            continue
        end
        
        if v1 isa Real
            @test v2 ≈ v1 nans=true
        elseif v1 isa AbstractVector
            @test all(isapprox.(v1, v2, nans=true, rtol=1e-10))
        end
    end
end


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
    @test_throws DomainError LilGuys.SurfaceDensityProfile([-1, 2, 3])
    @test_throws ArgumentError LilGuys.SurfaceDensityProfile([1])

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
        prof = LilGuys.SurfaceDensityProfile(radii, bins=bins)
        @test LilGuys.middle.(prof.log_Sigma) ≈ log_Sigma nans=true
        e = LilGuys.lower_error.(prof.log_Sigma)
        @test e[1:3] ≈ log_Sigma_em[1:3] atol=1e-8

        prof = LilGuys.SurfaceDensityProfile(radii, bins=bins, normalization=:mass)
        Mtot = 6
        e = LilGuys.lower_error.(prof.log_Sigma)
        @test LilGuys.middle.(prof.log_Sigma) ≈ log_Sigma .- log10(Mtot) nans=true atol=1e-8
        @test e[1:3] ≈ log_Sigma_em[1:3] atol=1e-8


        prof = LilGuys.SurfaceDensityProfile(radii, bins=bins, normalization=:central, bins_centre=1)
        Mtot = 2 / (π * 3.5^2)
        @test LilGuys.middle.(prof.log_Sigma) ≈ log_Sigma .- log10(Mtot) nans=true atol=1e-8
        e = LilGuys.lower_error.(prof.log_Sigma)
        @test e[1:3] ≈ log_Sigma_em[1:3] atol=1e-8


        @test_throws Exception LilGuys.SurfaceDensityProfile(radii, bins=bins, normalization=:jaberwocky)
    end


end


@testset "surface density of snapshot" begin
    @testset "simple test" begin
        x = [0., 0.5, 1.0, -1.0]
        y = [1.0, 0.0, -2.5, 3.0]
        z = [0.0, 1.0, 0.5, -1.0]


        positions = [x y z]'
        velocities = zeros(3, 4)
        masses = [1, 0, 1, 0]

        snap = LilGuys.Snapshot(positions, velocities, ones(4))
        snap.weights = masses

        prof = LilGuys.SurfaceDensityProfile(snap)

        prof = LilGuys.SurfaceDensityProfile(snap, x_vec = [0, 1, 0], y_vec = [0, 0, 1])

        prof = LilGuys.SurfaceDensityProfile(snap, R_units="arcmin")

        prof = LilGuys.SurfaceDensityProfile(snap, normalization=:none)
        @test false broken=true
    end
end


@testset "integration with exp profile unweighted" begin
    Σ(R) = exp(-R)/2π
    N = 10_000
    Rs = LilGuys.sample_surface_density(Σ, N, log_R=LinRange(-5, 5, 1000))

    obs = LilGuys.SurfaceDensityProfile(Rs, normalization=:none)
    
    Σs = 10 .^ LilGuys.middle.(obs.log_Sigma)
    mass_per_annulus = sum(Σs .* diff(π * 10 .^ 2obs.log_R_bins))

    @test sum(mass_per_annulus) ≈ N
    @test issorted(obs.log_R)
    @test sum(obs.counts) ≈ N

    R = 10 .^ obs.log_R
    sigma_exp = N * Σ.(R)
    log_sigma_exp = log10.(sigma_exp)
    
    @test_χ2 obs.log_Sigma log_sigma_exp

    Gamma_exp = -R
    err = maximum.(LilGuys.error_interval.(obs.Gamma))
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

        obs2 = LilGuys.SurfaceDensityProfile(filename)

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

    @testset "snapshot" begin
        theta = 2π*rand(N)
        x = @. Rs * cos(theta)
        y = @. Rs * sin(theta)
        z = Rs .* randn(N)
        positions = [x y z]'
        snap = Snapshot(positions, zeros(size(positions)), 1)
        obs2 = SurfaceDensityProfile(snap, normalization=:none, bins=obs.log_R_bins)
        @test obs.log_R ≈ obs2.log_R 
        @test obs.log_Sigma ≈ obs2.log_Sigma atol=0.08 nans=true

    end
end



@testset "integration (weighted)" begin
    N = 100_000
    radii = sqrt.(rand(N)) # uniform distribution
    bins = log10.(sqrt.(LinRange(minimum(radii), 1., 20)))

    obs = lguys.SurfaceDensityProfile(radii, bins=bins)
    sigma_exp = fill(N / π, length(obs.log_R))
    @test_χ2 obs.log_Sigma log10.(sigma_exp)

    prof = lguys.Plummer(r_s=0.2)
    weights = lguys.surface_density.(prof, radii)

    obs = lguys.SurfaceDensityProfile(radii, weights=weights, bins=bins)

    R = lguys.radii(obs)
    sigma_exp = N /π * surface_density.(prof, R)
    log_sigma_exp = log10.(sigma_exp)
    
    log_Sigma = lguys.log_surface_density(obs) .|> float
    log_Sigma_err = lguys.log_surface_density_err(obs) 

    # innermost bin severely underestimates uncertainty,
    # likely due to density gradient
    @test_χ2 log_Sigma[2:end] log_Sigma_err[2:end] log_sigma_exp[2:end]
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



@testset "integration with exp profile" begin
    Σ(r) = exp(-r)/2π
    N = 10_000
    r = LilGuys.sample_surface_density(Σ, N, log_R=LinRange(-5, 5, 1000))

    mass = 0.5 .+ 0.5rand(N)
    M = sum(mass)

    obs = LilGuys.SurfaceDensityProfile(r, normalization=:none, weights=mass,
        annotations = Dict("sigma" => π/3, "note" => "hi")
    )

    mass_per_annulus = 10 .^ obs.log_Sigma .* diff(π * 10 .^ 2obs.log_R_bins)
    @test sum(LilGuys.middle.(mass_per_annulus)) ≈ sum(mass)

    @test obs.annotations["sigma"] ≈ π/3
    @test obs.annotations["note"] == "hi"

    @test issorted(obs.log_R)
    @test sum(obs.counts) ≈ N

    r = 10 .^ obs.log_R
    sigma_exp = M * Σ.(r)

    err = maximum.(LilGuys.error_interval.(obs.log_Sigma))
    @test_χ2 LilGuys.middle.(obs.log_Sigma) err log10.(sigma_exp)

    Gamma_exp = -r
    err = maximum.(LilGuys.error_interval.(obs.Gamma))
    @test_χ2 LilGuys.middle.(obs.Gamma) err  Gamma_exp



    @testset "read/write" begin
        dir = mktempdir()

        filename = joinpath(dir, "test.toml")
        open(filename, "w") do f
            LilGuys.print(f, obs)
        end

        obs2 = LilGuys.SurfaceDensityProfile(filename)

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



@testset "scale" begin
    N = 100

    r = 10 .^ (0 .+ 0.4randn(N))
    m = rand(N)

    bins = LinRange(-1, 1, 5)
    prof_1 = LilGuys.SurfaceDensityProfile(r, weights=m, bins=bins, normalization=:none)

    M_scale = 0.232
    r_scale = 1.992
    prof_1 = LilGuys.scale(prof_1, r_scale, M_scale)

    prof_2 = LilGuys.SurfaceDensityProfile(r * r_scale, weights=m * M_scale, bins=bins .+ log10(r_scale), normalization=:none)

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



@testset "find longest consecutive_finite" begin
    @testset "simple cases" begin
        x = [1, NaN, 2, 3, 0, -1, Inf, 2, 3, NaN, 7, 8, 9]
        @test LilGuys.find_longest_consecutive_finite(x) == 3:6

        x = [0.23, -99.42, 5]
        @test LilGuys.find_longest_consecutive_finite(x) == 1:3

        @test LilGuys.find_longest_consecutive_finite([]) === 1:0
    end

end


@testset "edge_from_midpoint_filter" begin
    mids = 1:1
    @test LilGuys.edge_from_midpoint_filter(mids) == 1:2


    mids = 1:0
    @test LilGuys.edge_from_midpoint_filter(mids) == 1:1


    mids = 6:8
    @test LilGuys.edge_from_midpoint_filter(mids) == 6:9
end


@testset "filter_by_bin" begin
    prof = LilGuys.SurfaceDensityProfile(
        R_units="",
        log_R = [1.0, 1.5, 2.0],
        log_R_bins = [0.75, 1.25, 1.75, 2.25],
        counts = [1, 3, 2],
        log_Sigma = [Inf, 0.6, -0.2],
    )

    prof2 = LilGuys.filter_by_bin(prof, 1:1)

    @test prof2.log_R_bins ≈ [0.75, 1.25]

    @test prof2.log_R ≈ [1.0]
    @test prof2.counts ≈ [1]
    @test LilGuys.middle.(prof2.log_Sigma) ≈ [Inf]


    prof2 = LilGuys.filter_by_bin(prof, 2:3)

    @test prof2.log_R_bins ≈ [1.25, 1.75, 2.25]

    @test prof2.log_R ≈ [1.5, 2.0]
    @test prof2.counts ≈ [3, 2]
    @test LilGuys.middle.(prof2.log_Sigma) ≈ [0.6, -0.2]
end


@testset "Γ_max" begin
    @test false broken=true
end

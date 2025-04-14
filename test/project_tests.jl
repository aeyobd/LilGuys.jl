
@testset "to_sky" begin
    @testset "inverse" begin
        N = 100
        snap = lguys.Snapshot(positions=100randn(3, N), velocities=1randn(3, N), masses=ones(N), index=1:N, header=lguys.make_default_header(1, N), potential=-ones(N))

        sky = lguys.to_sky(snap)
        gc = lguys.transform.(lguys.Galactocentric, sky)

        positions = [[g.x, g.y, g.z] for g in gc]
        velocities = [[g.v_x, g.v_y, g.v_z] for g in gc]
        positions = hcat(positions...)
        velocities = hcat(velocities...)

        velocities ./= V2KMS

        @test positions ≈ snap.positions
        @test velocities ≈ snap.velocities
    end
end


@testset "to_gaia" begin
    N = 100
    pos = 100randn(3, N)
    vel = 1randn(3, N)
    masses = ones(N)
    weights = rand(N)

    snap = lguys.Snapshot(positions=pos, velocities=vel, masses=masses, weights=weights,
        index=1:N, header=lguys.make_default_header(1, N), potential=-ones(N))

    gaia = lguys.to_gaia(snap, add_centre=false)

    @testset "integration" begin
        obs = lguys.to_sky(snap) 

        @test size(gaia, 1) == N 
        @test gaia.weights ≈ weights
        @test ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"] ⊆ names(gaia)

        @test gaia.distance ≈ [o.distance for o in obs]
        @test gaia.ra ≈ [o.ra for o in obs]

        # test chi2, r_ell, etc
        @test false broken=true
    end

    @testset "add centre" begin
        # test add centre
        gaia_cen = lguys.to_gaia(snap, add_centre=true)
        @test size(gaia_cen, 1) == N + 1
        @test gaia_cen.weights[1] == 0
        @test gaia_cen.distance[1] == 0 broken=true

    end
    
    @testset "filters" begin
        @test false broken=true 
    end
end

@testset "to_frame" begin
    @test false broken=true 
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




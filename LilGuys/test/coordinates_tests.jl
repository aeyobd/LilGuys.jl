

@testset "cartesian initialization" begin

    position = rand(3)
    velocity = rand(3)

    x, y, z = position
    v_x, v_y, v_z = velocity

    struct TestFrame <: lguys.CoordinateFrame end

    p1 = lguys.Cartesian{TestFrame}(position, velocity)
    p2  = lguys.Cartesian{TestFrame}(x=x, y=y, z=z, v_x=v_x, v_y=v_y, v_z=v_z)

    @test fieldnames(lguys.Cartesian{TestFrame}) == (:x, :y, :z, :v_x, :v_y, :v_z, :coord)

    for sym in fieldnames(lguys.Cartesian{TestFrame})
        @test getfield(p1, sym) == getfield(p2, sym)
    end

    @test lguys.position_of(p1) == position
    @test lguys.velocity_of(p1) == velocity

    @test lguys.position_of(p2) == position
    @test lguys.velocity_of(p2) == velocity
end


@testset "phase point repr" begin
    position = [0.23523, 1.2342, π]
    velocity = [exp(1), 0, -234.234]

    p = lguys.Cartesian{TestFrame}(position, velocity)

    @test repr(p) == "TestFrame point at (0.24, 1.23, 3.14) kpc, (2.72, 0.00, -234.23) km/s"
end


@testset "sky coord repr" begin
    ra = 0.23523
    dec = -1.2342

    sc = lguys.ICRS(ra=ra, dec=dec)

    @test repr(sc) == "ICRS at (0.24, -1.23) deg"
end



@testset "initialize sky coordinates" begin
    frames = [lguys.ICRS, lguys.GSR]

    for frame in frames
        ra = rand() * 360
        dec = rand() * 180 - 90

        distance = rand() * 1000
        pmra = rand() * 100
        pmdec = rand() * 100
        radial_velocity = rand() * 100

        sc = frame(ra=ra, dec=dec, distance=distance, pmra=pmra, pmdec=pmdec, radial_velocity=radial_velocity)

        @test sc.ra ≈ ra
        @test sc.dec ≈ dec
        @test sc.distance ≈ distance
        @test sc.pmra ≈ pmra
        @test sc.pmdec ≈ pmdec
        @test sc.radial_velocity ≈ radial_velocity
    end
end


@testset "initialize cartesian coordinates" begin
    frames = [lguys.Cartesian{lguys.ICRS}, lguys.Cartesian{lguys.GSR}, lguys.Galactocentric]

    for frame in frames
        position = rand(3)
        velocity = rand(3)

        x, y, z = position
        v_x, v_y, v_z = velocity

        p = frame(position, velocity)

        @test p.x ≈ x
        @test p.y ≈ y
        @test p.z ≈ z
        @test p.v_x ≈ v_x
        @test p.v_y ≈ v_y
        @test p.v_z ≈ v_z
    end
end

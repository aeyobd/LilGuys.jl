@testset "Point3D" begin 
    x, y, z = 1, 2, 3
    p = lguys.Point3D(x, y, z)

    @test p.x == x
    @test p.y == y
    @test p.z == z

    @testset "promotion" begin
        p1 = lguys.Point3D(25, 0.2, π)
        @test typeof(p1) == lguys.Point3D{Float64}
        @test p1.x == 25.0
        @test p1.y == 0.2
        @test p1.z ≈ π
    end
end


@testset "Point6D" begin
    p = lguys.Point6D(1, 2, 3, 4, 5, 6)
    @test p.x == 1
    @test p.y == 2
    @test p.z == 3
    @test p.v_x == 4
    @test p.v_y == 5
    @test p.v_z == 6

    @testset "promotion" begin
        p1 = lguys.Point6D(25, 0.2, π, 1//2, 10, Inf)
        @test typeof(p1) == lguys.Point6D{Float64}
        @test p1.x == 25.0
        @test p1.y == 0.2
        @test p1.z ≈ π
        @test p1.v_x == 0.5
        @test p1.v_y == 10.0
        @test p1.v_z === Inf
    end
end



@testset "cartesian initialization" begin

    position = rand(3)
    velocity = rand(3)

    x, y, z = position
    v_x, v_y, v_z = velocity

    struct TestFrame <: lguys.CoordinateFrame end

    p1 = lguys.Cartesian{TestFrame}(position, velocity)
    p2  = lguys.Cartesian{TestFrame, Float64}(x, y, z, v_x, v_y, v_z)

    # @test fieldnames(lguys.Cartesian{TestFrame, Float64}) == (:x, :y, :z, :v_x, :v_y, :v_z, :coord)

    for sym in fieldnames(lguys.Cartesian{TestFrame, Float64})
        @test getfield(p1, sym) == getfield(p2, sym)
    end

    @test lguys.position(p1) == position
    @test lguys.velocity(p1) == velocity

    @test lguys.position(p2) == position
    @test lguys.velocity(p2) == velocity
end


@testset "phase point repr" begin
    position = [0.23523, 1.2342, π]
    velocity = [exp(1), 0, -234.234]

    p = lguys.Cartesian{TestFrame}(position, velocity)

    @test repr(p) == "TestFrame point at (0.24, 1.23, 3.14) kpc, (2.72, 0.00, -234.23) km/s" broken=true
end


@testset "sky coord repr" begin
    ra = 0.23523
    dec = -1.2342

    sc = lguys.ICRS(ra=ra, dec=dec)

    @test repr(sc) == "ICRS{Float64} at (0.24, -1.23) deg" broken=true
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
    frames = [lguys.Cartesian{lguys.ICRS, Float64}, lguys.Cartesian{lguys.GSR, Float64}, lguys.Galactocentric]

    for frame in frames
        position = rand(3)
        velocity = rand(3)

        x, y, z = position
        v_x, v_y, v_z = velocity

        p = frame(x, y, z, v_x, v_y, v_z)

        @test p.x ≈ x
        @test p.y ≈ y
        @test p.z ≈ z
        @test p.v_x ≈ v_x
        @test p.v_y ≈ v_y
        @test p.v_z ≈ v_z
    end
end


@testset "rand_coords" begin
    obs = lguys.ICRS(ra=15, dec=-33, distance=80, pmra=-0.1, pmdec=0.2, radial_velocity=50)
    err = lguys.ICRS(ra=0.1, dec=0.1, distance=1, pmra=0.01, pmdec=0.01, radial_velocity=0.4)

    N = 1000

    coords = lguys.rand_coords(obs, err, N)

    @test length(coords) == N
    @test lguys.mean([c.ra for c in coords]) ≈ obs.ra atol=0.1
    @test lguys.mean([c.dec for c in coords]) ≈ obs.dec atol=0.1
    @test lguys.mean([c.distance for c in coords]) ≈ obs.distance atol=1
    @test lguys.mean([c.pmra for c in coords]) ≈ obs.pmra atol=0.01
    @test lguys.mean([c.pmdec for c in coords]) ≈ obs.pmdec atol=0.01
    @test lguys.mean([c.radial_velocity for c in coords]) ≈ obs.radial_velocity atol=0.4

    @test lguys.std([c.ra for c in coords]) ≈ err.ra atol=0.1
    @test lguys.std([c.dec for c in coords]) ≈ err.dec atol=0.1
    @test lguys.std([c.distance for c in coords]) ≈ err.distance atol=1
    @test lguys.std([c.pmra for c in coords]) ≈ err.pmra atol=0.01
    @test lguys.std([c.pmdec for c in coords]) ≈ err.pmdec atol=0.01
    @test lguys.std([c.radial_velocity for c in coords]) ≈ err.radial_velocity atol=0.4
end


@testset "rand_split_gaussian" begin

end


@testset "rand_coords" begin
    N = 1000
    df = Dict(
        "ra" => 15,
        "ra_err" => 0.1,
        "dec" => -33,
        "dec_err" => 0.1,
        "distance_modulus" => 18,
        "distance_modulus_err" => 0.3,
        "pmra" => -0.1,
        "pmra_err" => 0.01,
        "pmdec" => 0.2,
        "pmdec_err" => 0.01,
        "radial_velocity" => 50,
        "radial_velocity_err" => 0.4
       )

    dist = lguys.dm2kpc(df["distance_modulus"])
    dist_err = dist * df["distance_modulus_err"] / 5 * log(10)
    coords = lguys.rand_coords(df, N)

    @test length(coords) == N
    @test lguys.mean([c.ra for c in coords]) ≈ df["ra"] atol=0.1
    @test lguys.mean([c.dec for c in coords]) ≈ df["dec"] atol=0.1
    @test lguys.mean([c.distance for c in coords]) ≈ dist atol=1
    @test lguys.mean([c.pmra for c in coords]) ≈ df["pmra"] atol=0.01
    @test lguys.mean([c.pmdec for c in coords]) ≈ df["pmdec"] atol=0.01
    @test lguys.mean([c.radial_velocity for c in coords]) ≈ df["radial_velocity"] atol=0.4

    @test lguys.std([c.ra for c in coords]) ≈ df["ra_err"] atol=0.1
    @test lguys.std([c.dec for c in coords]) ≈ df["dec_err"] atol=0.1
    @test lguys.std([c.distance for c in coords]) ≈ dist_err atol=1
    @test lguys.std([c.pmra for c in coords]) ≈ df["pmra_err"] atol=0.01
    @test lguys.std([c.pmdec for c in coords]) ≈ df["pmdec_err"] atol=0.01
    @test lguys.std([c.radial_velocity for c in coords]) ≈ df["radial_velocity_err"] atol=0.4
end



@testset "coords from file" begin
end

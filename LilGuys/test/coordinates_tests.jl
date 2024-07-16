

@testset "phase point initialization" begin

    position = rand(3)
    velocity = rand(3)

    x, y, z = position
    v_x, v_y, v_z = velocity

    p1 = lguys.PhasePoint{Test}(position, velocity)
    p2  = lguys.PhasePoint{Test}(x, y, z, v_x, v_y, v_z)

    @test fieldnames(lguys.PhasePoint{Test}) == (:x, :y, :z, :v_x, :v_y, :v_z)

    for sym in fieldnames(lguys.PhasePoint{Test})
        @test getfield(p1, sym) == getfield(p2, sym)
    end

    @test lguys.get_position(p1) == position
    @test lguys.get_velocity(p1) == velocity

    @test lguys.get_position(p2) == position
    @test lguys.get_velocity(p2) == velocity
end


@testset "phase point repr" begin
    position = [0.23523, 1.2342, π]
    velocity = [exp(1), 0, -234.234]

    p = lguys.PhasePoint{Test}(position, velocity)

    @test repr(p) == "Test point at (0.24, 1.23, 3.14) kpc, (2.72, 0.00, -234.23) km/s"
end


@testset "sky coord repr" begin
    ra = 0.23523
    dec = -1.2342

    sc = lguys.SkyCoord{Test}(ra=ra, dec=dec)

    @test repr(sc) == "Test at (0.24, -1.23) deg"
end



@testset "initialize coordinates" begin
    frames = [lguys.ICRS, lguys.ICRS_Cartesian, lguys.GSR, lguys.GSR_Cartesian, lguys.Galactocentric]

    for frame in frames
    end
end

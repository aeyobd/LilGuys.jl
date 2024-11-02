sun_gc = [-8.121973366, 0, 0.0208] # kpc, uses d cos(θ), 0, d sin(θ) where θ  = arcsin(z_sun / d_sun)

coordinate_systems = [
    lguys.ICRS, lguys.Cartesian{lguys.ICRS},
    lguys.GSR, lguys.Cartesian{lguys.GSR},
    lguys.Galactocentric 
   ]


@testset "sky to cartesian" begin
    icrs = lguys.ICRS(ra=0., dec=0, distance=1,
                      pmra=0, pmdec=0, radial_velocity=0)

    cart = lguys.transform(lguys.Cartesian{lguys.ICRS}, icrs)

    @test cart.x ≈ 1 rtol=1e-2
    @test cart.y ≈ 0 rtol=1e-2
    @test cart.z ≈ 0 rtol=1e-2
    @test cart.v_x ≈ 0 rtol=1e-2
    @test cart.v_y ≈ 0 rtol=1e-2
    @test cart.v_z ≈ 0 rtol=1e-2

end

@testset "cartesian to sky" begin

    cart = lguys.Cartesian{lguys.ICRS}(x=1., y=0, z=0,
                                      v_x=0, v_y=0, v_z=0)
    icrs = lguys.transform(lguys.ICRS, cart)

    @test icrs.ra ≈ 0 rtol=1e-2
    @test icrs.dec ≈ 0 rtol=1e-2
    @test icrs.distance ≈ 1 rtol=1e-2
    @test icrs.pmra ≈ 0 rtol=1e-2
    @test icrs.pmdec ≈ 0 rtol=1e-2
    @test icrs.radial_velocity ≈ 0 rtol=1e-2
end


@testset "helio to galcen: Sag A*" begin
    gc = lguys.ICRS(ra = 266.4051, dec=-28.936175, distance=8.122,
                             pmra=-3.151, pmdec=-5.547, radial_velocity=-12.9)

    phase = lguys.transform(lguys.Galactocentric, gc)

    @test lguys.position_of(phase) ≈ [0,0,0] atol=1e-2
    @test phase.v_x ≈ 0 atol=0.1
    @test phase.v_y ≈ 0 atol=0.1
    @test phase.v_z ≈ 0 atol=0.1



    sun = lguys.ICRS(ra = 0., dec=-0, distance=0,
                             pmra=0, pmdec=0, radial_velocity=0)
    phase = lguys.transform(lguys.Galactocentric, sun)
    @test lguys.position_of(phase) ≈ sun_gc rtol=3e-3
    @test phase.v_x ≈12.9 atol=0.1
    @test phase.v_y ≈ 245.6 atol=0.1
    @test phase.v_z ≈ 7.78 atol=0.1
end


@testset "galcen to helio: Sag A*" begin 
    g = lguys.Galactocentric(x=0., y=0, z=0,
                             v_x=0, v_y=0, v_z=0)
    obs = lguys.transform(lguys.ICRS, g)

    @test obs.ra ≈ 266.4168166 rtol=1e-2
    @test obs.dec ≈ -29.00782 rtol=1e-2
    @test obs.distance ≈ 8.122 rtol=1e-2
    @test obs.pmra ≈ -3.151 rtol=1e-2
    @test obs.pmdec ≈ -5.547 rtol=1e-2
    @test obs.radial_velocity ≈ -12.9 atol=0.1
end




@testset "helio to galcen: inverse" begin
    N = 100
    obs = [lguys.ICRS(ra=360rand(), dec=-90 + 180rand(), distance=2*rand(),
                                pmra=10*randn(), pmdec=10*randn(), radial_velocity=100*randn())
             for _ in 1:N]

    obs2 = lguys.transform.(lguys.ICRS, lguys.transform.(lguys.Galactocentric, obs))

    for i in 1:N
        p = obs[i]
        q = obs2[i]
        @test p.ra ≈ q.ra rtol=1e-2
        @test p.dec ≈ q.dec rtol=1e-2
        @test p.distance ≈ q.distance rtol=1e-2
        @test p.pmra ≈ q.pmra rtol=1e-2
        @test p.pmdec ≈ q.pmdec rtol=1e-2
        @test p.radial_velocity ≈ q.radial_velocity rtol=1e-2
    end

end


@testset "GSR to galcen: Sag A*" begin
    gc = lguys.GSR(ra = 266.4051, dec=-28.936175, distance=8.122,
                   pmra=0., pmdec=0., radial_velocity=0.)

    phase = lguys.transform(lguys.Galactocentric, gc)

    @test lguys.position_of(phase) ≈ [0,0,0] atol=1e-2
    @test lguys.velocity_of(phase) ≈ [0,0,0] atol=1e-2



    static_sun = lguys.GSR(ra=0., dec=0, distance=0,
                             pmra=0, pmdec=0, radial_velocity=0)
    phase = lguys.transform(lguys.Galactocentric, static_sun)

    frame = phase.frame
    theta = asin(frame.z_sun / frame.d)
    @test lguys.position_of(phase) ≈ [-8.122*cos(theta), 0, 8.122 * sin(theta)] atol=3e-3
    @test lguys.velocity_of(phase) ≈ [0, 0, 0] atol=3e-3


    gsr = lguys.GSR(ra=266.4051, dec=-28.936175, distance=8.122,
                    pmra=-3.151, pmdec=-5.547, radial_velocity=-12.9)

    gc = lguys.transform(lguys.Galactocentric, gsr)
    @test lguys.position_of(gc) ≈ [0,0,0] atol=1e-2
    @test lguys.velocity_of(gc) ≈ - frame.v_sun atol = 0.2
end


@testset "galcen to GSR: Sag A*" begin 
    g = lguys.Galactocentric(zeros(3), zeros(3))
    obs = lguys.transform(lguys.GSR, g)

    @test obs.ra ≈ 266.4168166 rtol=1e-2
    @test obs.dec ≈ -29.00782 rtol=1e-2
    @test obs.distance ≈ 8.122 rtol=1e-2
    @test obs.pmra ≈ 0 atol=1e-2
    @test obs.pmdec ≈ 0 atol=1e-2
    @test obs.radial_velocity ≈ 0 atol=0.1
end


@testset "ICRS to GSR: example" begin
    # from astropy
    icrs = lguys.ICRS(ra=258.58356362, dec=14.55255619, distance=10.0,
                      pmra=0, pmdec=0, radial_velocity=-16.1)

    gsr = lguys.transform(lguys.GSR, icrs)

    @test gsr.ra ≈ icrs.ra rtol=1e-2
    @test gsr.dec ≈ icrs.dec rtol=1e-2

    @test gsr.distance ≈ icrs.distance rtol=1e-2
    @test gsr.radial_velocity ≈ 123.30460087379765 rtol=1e-2
end


function inverse_test(frame1::Type{<:lguys.AbstractCartesian}, frame2)
    N = 10

    p1 = [frame1(x=rand(), y=rand(), z=rand(),
                 v_x=rand(), v_y=rand(), v_z=rand())
             for _ in 1:N]

    p2 = lguys.transform.(frame1, lguys.transform.(frame2, p1))

    for field in [:x, :y, :z, :v_x, :v_y, :v_z]
        @test getproperty.(p1, field) ≈ getproperty.(p2, field) rtol=1e-2
    end

end

function inverse_test(frame1::Type{<:lguys.AbstractSkyCoord}, frame2)
    N = 10

    p1 = [frame1(ra=360rand(), dec=-90 + 180rand(), distance=100*rand(),
                 pmra=10*randn(), pmdec=10*randn(), radial_velocity=100*randn())
             for _ in 1:N]

    p2 = lguys.transform.(frame1, lguys.transform.(frame2, p1))

    for field in [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]
        @test getproperty.(p1, field) ≈ getproperty.(p2, field) rtol=1e-2
    end

end


@testset "inverse" begin 
    for frame1 in coordinate_systems
        for frame2 in coordinate_systems
            @testset "inverse - $frame1 -> $frame2 -> $frame1" begin
                inverse_test(frame1, frame2)
            end

        end
    end
end


@testset "galactocentric frame settings" begin
    @test false broken=true
end

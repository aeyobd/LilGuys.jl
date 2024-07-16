coordinate_systems = [
    lguys.ICRS, lguys.ICRS_Cartesian,
    lguys.GSR, lguys.GSR_Cartesian,
    lguys.Galactocentric 
   ]

@testset "helio to galcen: Sag A*" begin
    gc = lguys.ICRS(ra = 266.4051, dec=-28.936175, distance=8.122,
                             pmra=-3.151, pmdec=-5.547, radial_velocity=-12.9)

    phase = lguys.transform(lguys.Galactocentric, gc)

    @test lguys.get_position(phase) ≈ [0,0,0] atol=1e-2
    @test phase.v_x ≈ 0 atol=0.1
    @test phase.v_y ≈ 0 atol=0.1
    @test phase.v_z ≈ 0 atol=0.1



    sun = lguys.ICRS(ra = 0, dec=-0, distance=0,
                             pmra=0, pmdec=0, radial_velocity=0)
    phase = lguys.transform(lguys.Galactocentric, sun)
    @test lguys.get_position(phase) ≈ [-8.122, 0, 0] rtol=3e-3
    @test phase.v_x ≈12.9 atol=0.1
    @test phase.v_y ≈ 245.6 atol=0.1
    @test phase.v_z ≈ 7.78 atol=0.1
end


@testset "galcen to helio: Sag A*" begin 
    g = lguys.Galactocentric(zeros(3), zeros(3))
    obs = lguys.transform(lguys.ICRS, g)

    @test obs.ra ≈ 266.4168166 rtol=1e-2
    @test obs.dec ≈ -29.00782 rtol=1e-2
    @test obs.distance ≈ 8.122 rtol=1e-2
    @test obs.pmra ≈ -3.151 rtol=1e-2
    @test obs.pmdec ≈ -5.547 rtol=1e-2
    @test obs.radial_velocity ≈ -12.9 atol=0.1
end


# @testset "inverse" begin 
#     for frame1 in coordinate_systems
#         for frame2 in coordinate_systems
#             N = 100
#             phase = [frame1(
#                      for _ in 1:N]
#             phase2 = lguys.transform.(lguys.Galactocentric, lguys.transform.(lguys.ICRS, phase))
# 
#             for i in 1:N
#                 p = phase[i]
#                 q = phase2[i]
#                 for field in fieldnames(frame1)
#                     @test p.(field) ≈ q.(field) rtol=1e-2
#                 end
#             end
#         end
#     end
# end


@testset "helio to galcen: inverse" begin
    N = 100
    obs = [lguys.ICRS(360rand(), -90 + 180rand(), 2*rand(),
                                10*randn(), 10*randn(), 10*randn())
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
                   pmra=0, pmdec=0, radial_velocity=0)

    phase = lguys.transform(lguys.Galactocentric, gc)

    @test lguys.get_position(phase) ≈ [0,0,0] atol=1e-2
    @test lguys.get_velocity(phase) ≈ [0,0,0] atol=1e-2



    static_sun = lguys.GSR(ra = 0, dec=0, distance=0,
                             pmra=0, pmdec=0, radial_velocity=0)
    phase = lguys.transform(lguys.Galactocentric, static_sun)

    frame = phase.frame
    theta = asin(frame.z_sun / frame.d)
    @test lguys.get_position(phase) ≈ [-8.122*cos(theta), 0, 8.122 * sin(theta)] atol=3e-3
    @test lguys.get_velocity(phase) ≈ [0, 0, 0] atol=3e-3


    gsr = lguys.GSR(ra=266.4051, dec=-28.936175, distance=8.122,
                    pmra=-3.151, pmdec=-5.547, radial_velocity=-12.9)

    gc = lguys.transform(lguys.Galactocentric, gsr)
    @test lguys.get_position(gc) ≈ [0,0,0] atol=1e-2
    @test lguys.get_velocity(gc) ≈ - frame.v_sun atol = 0.2
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

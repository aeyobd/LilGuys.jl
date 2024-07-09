
@testset "helio to galcen: Sag A*" begin
    gc = lguys.ICRS(ra = 266.4051, dec=-28.936175, distance=8.122,
                             pmra=-3.151, pmdec=-5.547, radial_velocity=-12.9)

    phase = lguys.transform(lguys.Galactocentric, gc)

    @test all(abs.(phase.position) .< 1e-2)
    @test phase.velocity ≈ [0,0,0] atol=0.2 # TODO this is really high



    sun = lguys.ICRS(ra = 0, dec=-0, distance=0,
                             pmra=0, pmdec=0, radial_velocity=0)
    phase = lguys.transform(lguys.Galactocentric, sun)
    @test phase.position ≈ [-8.122, 0, 0] rtol=3e-3
    @test phase.velocity ≈ [12.9, 245.6, 7.78] rtol=3e-3
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


@testset "galcen to helio: inverse" begin
    N = 100
    phase = [lguys.Galactocentric( 20*randn(3), 100*randn(3))
             for _ in 1:N]
    phase2 = lguys.transform.(lguys.Galactocentric, lguys.transform.(lguys.ICRS, phase))

    for i in 1:N
        p = phase[i]
        q = phase2[i]
        @test p.position ≈ q.position rtol=1e-2
        @test p.velocity ≈ q.velocity rtol=1e-2
    end

end


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

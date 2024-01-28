
@testset "test to galcen" begin
    gc = lguys.Observation(ra = 266.4051, dec=-28.936175, distance=8.122,
                             pm_ra=-3.151, pm_dec=-5.547, radial_velocity=-12.9)

    phase = lguys.transform(lguys.Galactocentric, gc)

    @test all(abs.(phase.position) .< 1e-2)
    @test phase.velocity ≈ [0,0,0] atol=0.2 # TODO this is really high



    sun = lguys.Observation(ra = 0, dec=-0, distance=0,
                             pm_ra=0, pm_dec=0, radial_velocity=0)
    phase = lguys.transform(lguys.Galactocentric, sun)
    @test phase.position ≈ [-8.122, 0, 0] rtol=3e-3
    @test phase.velocity ≈ [12.9, 245.6, 7.78] rtol=3e-3
end


@testset "test to geocen" begin 
    g = lguys.Galactocentric(zeros(3), zeros(3))
    obs = lguys.transform(lguys.Observation, g)

    @test obs.ra ≈ 266.4168166 rtol=1e-2
    @test obs.dec ≈ -29.00782 rtol=1e-2
    @test obs.distance ≈ 8.122 rtol=1e-2
    @test obs.pm_ra ≈ -3.151 rtol=1e-2
    @test obs.pm_dec ≈ -5.547 rtol=1e-2
    @test obs.radial_velocity ≈ -12.9 atol=0.1
end


@testset "test inverse" begin
    N = 100
    phase = [lguys.Galactocentric( 20*randn(3), 100*randn(3))
             for _ in 1:N]
    phase2 = lguys.transform.(lguys.Galactocentric, lguys.transform.(lguys.Observation, phase))

    for i in 1:N
        p = phase[i]
        q = phase2[i]
        @test p.position ≈ q.position rtol=1e-2
        @test p.velocity ≈ q.velocity rtol=1e-2
    end

end


@testset "test inverse 2" begin
    N = 100
    obs = [lguys.Observation(360rand(), -90 + 180rand(), 2*rand(),
                                10*randn(), 10*randn(), 10*randn())
             for _ in 1:N]

    obs2 = lguys.transform.(lguys.Observation, lguys.transform.(lguys.Galactocentric, obs))

    for i in 1:N
        p = obs[i]
        q = obs2[i]
        @test p.ra ≈ q.ra rtol=1e-2
        @test p.dec ≈ q.dec rtol=1e-2
        @test p.distance ≈ q.distance rtol=1e-2
        @test p.pm_ra ≈ q.pm_ra rtol=1e-2
        @test p.pm_dec ≈ q.pm_dec rtol=1e-2
        @test p.radial_velocity ≈ q.radial_velocity rtol=1e-2
    end

end

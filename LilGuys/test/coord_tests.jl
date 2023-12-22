import LilGuys


@testset "test to galcen" begin
    gc = LilGuys.Observation(ra = 266.4168166, dec=-29.00782, distance=8.29,
                             pm_ra=0, pm_dec=0, radial_velocity=0)

    phase = LilGuys.to_galcen(gc)

    @test all(abs.(phase.pos) .< 1e-2)
end


@testset "test to geocen" begin 
    g = LilGuys.PhasePoint(zeros(3), zeros(3))
    obs = LilGuys.to_sky(g)

    @test obs.ra ≈ 266.4168166 rtol=1e-2
    @test obs.dec ≈ -29.00782 rtol=1e-2
    @test obs.distance ≈ 8.29 rtol=1e-2

end

@testset "test inverse" begin
    N = 100
    phase = [LilGuys.PhasePoint( 20*randn(3), 100*randn(3))
             for _ in 1:N]
    phase2 = LilGuys.to_galcen.(LilGuys.to_sky.(phase))

    for i in 1:N
        p = phase[i]
        q = phase2[i]
        @test p.pos ≈ q.pos rtol=1e-2
        @test p.vel ≈ q.vel rtol=1e-2
    end

end

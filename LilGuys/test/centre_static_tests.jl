
@testset "centre of mass" begin

end

@testset "potential weighted" begin
end


@testset "centroid snap" begin
    N = 10
    pos = transpose(hcat(ones(N), zeros(N), -ones(N)))
    vel = transpose(hcat(zeros(N), 2*ones(N), zeros(N)))
    m = lguys.ConstVector(1., N)
    snap = lguys.Snapshot(pos, vel, m)

    cen = lguys.Centres.centre_of_mass(snap)

    @test cen.position ≈ [1., 0., -1.] 
    @test cen.velocity ≈ [0., 2., 0.] 
    @test cen.position_err ≈ 0 atol=1e-10
    @test cen.velocity_err ≈ 0 atol=1e-10
end


@testset "centroid snap weighted" begin
    N = 100000
    pos = randn(3, N)
    vel = randn(3, N)
    m = rand(N)
    snap = lguys.Snapshot(pos, vel, m)

    cen = lguys.Centres.weighted_centre(snap, m)

    @test cen.position ≈ [0., 0., 0.] atol=1e-2
    @test cen.velocity ≈ [0., 0., 0.] atol=1e-2
    @test cen.position_err ≈ 1/sqrt(N) atol=1e-2
    @test cen.velocity_err ≈ 1/sqrt(N) atol=1e-2
end


@testset "static centre potential" begin
    N = 100000
    pos = randn(3, N)
    vel = randn(3, N)
    snap = lguys.Snapshot(pos, vel, ones(N))
    snap.Φs = lguys.calc_radial_discrete_Φ(snap)

    cen1 = lguys.Centres.weighted_centre(snap, -snap.Φs)
    cen2 = lguys.Centres.centre_weighted_potential(snap)

    s2 = lguys.StaticState(snap; method="potential")
    cen3 = lguys.calc_centre!(s2, snap).centre
    
    @test cen1.position ≈ cen2.position atol=1e-2
    @test cen1.velocity ≈ cen2.velocity atol=1e-2
    @test cen1.position_err ≈ cen2.position_err atol=1e-2
    @test cen1.velocity_err ≈ cen2.velocity_err atol=1e-2
    
    @test cen1.position ≈ cen3.position atol=1e-2
    @test cen1.velocity ≈ cen3.velocity atol=1e-2
    @test cen1.position_err ≈ cen3.position_err atol=1e-2
    @test cen1.velocity_err ≈ cen3.velocity_err atol=1e-2
end

@testset "centre_potential_percen" begin
    N = 100000
    pos = randn(3, N)
    vel = randn(3, N)

    m = rand(N)
    snap2 = lguys.Snapshot(pos, vel, m)
    cen = lguys.Centres.centre_potential_percen(snap2)

    snap2.Φs = lguys.calc_radial_discrete_Φ(snap2)

    snap = snap2[snap2.Φs .< lguys.percentile(snap2.Φs, 5)]

    cen2 = lguys.Centres.centre_of_mass(snap)


    @test cen.position ≈ cen2.position atol=1e-2
    @test cen.velocity ≈ cen2.velocity atol=1e-2
    @test cen.position_err ≈ cen2.position_err atol=1e-2
    @test cen.velocity_err ≈ cen2.velocity_err atol=1e-2
end

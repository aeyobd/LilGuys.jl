@testset "shrinking spheres params: initialization" begin
    positions = randn(3, 100)
    x0 = lguys.centroid(positions)

    params = lguys.Centres._ShrinkingSpheresParams(positions)

    @test params.N == 100
    @test params.x0 == x0


    @test maximum(lguys.calc_r(positions, x0)) == params.r_cut_0

end


@testset "shrinking spheres params: N_min" begin
    N = 200
    positions = randn(3, N)

    N_min = 50
    f_min = 0.101
    params = lguys.Centres._ShrinkingSpheresParams(positions, N_min=N_min, f_min=f_min)

    @test params.N_min == 50


    N_min = 12
    f_min = 0.2001
    params = lguys.Centres._ShrinkingSpheresParams(positions, N_min=N_min, f_min=f_min)

    @test params.N_min == 40

    
    params = lguys.Centres._ShrinkingSpheresParams(positions, N_min=N_min)

    @test params.N_min == N_min


    f_min = 0.6001
    params = lguys.Centres._ShrinkingSpheresParams(positions, f_min=f_min)

    @test params.N_min == 120


    params = lguys.Centres._ShrinkingSpheresParams(positions)

    @test params.N_min == 100
end


@testset "shrinking spheres params validation" begin 
    params = lguys.Centres._ShrinkingSpheresParams(randn(3, 100))

end


@testset "shrinking spheres: one point" begin 
    positions = [2.;0.12;0.25;;]
    params = lguys.Centres._ShrinkingSpheresParams(positions)
    println(params.x0)
    cen, filt = lguys.shrinking_spheres(positions)
    @test cen ≈ [2.;0.12;0.25] rtol=1e-10
    @test sum(filt) == 1
end



@testset "shrinking spheres centre: trivial" begin
    positions = [2.;0.12;0.25;;]
    velocities = [0.1;0.1;0.4;;]
    masses = [1.]

    snap = lguys.Snapshot(positions, velocities, masses)
    snap.Φs = [0.]
    snap.accelerations = [0.;0.;0.;;]

    state = lguys.SS_State(snap, verbose=true)
    
    println(snap.positions)
    println(state)

    lguys.Centres.calc_centre!(state, snap)
    cen = state.centre

    @test cen.position ≈ positions[:, 1]
    @test cen.velocity ≈ velocities[:, 1]
end





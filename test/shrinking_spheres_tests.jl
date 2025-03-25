using Logging

@testset "shrinking spheres params: initialization" begin
    positions = randn(3, 100)
    x0 = lguys.centroid(positions)

    params = lguys.Centres._ShrinkingSpheresParams(positions)

    @test params.N == 100
    @test params.x0 == x0


    @test maximum(lguys.radii(positions, x0)) == params.r_cut_0

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
    cen, filt = lguys.shrinking_spheres(positions)
    @test cen ≈ [2.;0.12;0.25] rtol=1e-10
    @test sum(filt) == 1
end



@testset "shrinking spheres centre: trivial" begin
    positions = [2.;0.12;0.25;;]
    velocities = [0.1;0.1;0.4;;]
    masses = [1.]

    snap = lguys.Snapshot(positions, velocities, masses)
    snap.potential = [0.]
    snap.accelerations = [0.;0.;0.;;]

    state = lguys.SS_State(snap)
    

    lguys.Centres.calc_centre!(state, snap)
    cen = state.centre

    @test cen.position ≈ positions[:, 1]
    @test cen.velocity ≈ velocities[:, 1]
end





@testset "shrinking spheres NFW" begin
    # check converges for many different initial positions
    r_s=3.0
    M_s=0.8

    N = 10_000
    v_s = sqrt(M_s / r_s)
    nfw = LilGuys.TruncNFW(r_s=r_s,  M_s=M_s, trunc=10)
    snap = LilGuys.sample_potential(nfw, N)

    snap.potential = LilGuys.potential_spherical_discrete(snap)

    # TODO; these sem reasonable but not sure how to properly quantify this
    dx = r_s * 0.23 / sqrt(100)
    dv = v_s * 0.025 / sqrt(100)

    @testset  "zeroed" begin
        test_logger = TestLogger()

        state = lguys.SS_State(snap)

        with_logger(test_logger) do
            LilGuys.Centres.calc_centre!(state, snap)
        end
        cen = state.centre

        @test cen.position ≈ [0.,0.,0.] atol=0.5
        @test cen.velocity ≈ [0.,0.,0.] atol=0.005
        @test cen.position_err ≈ dx rtol=0.5
        @test cen.velocity_err ≈ dv rtol=0.5
    end


    @testset "shifted centre" begin
        state = lguys.SS_State(snap)
        state.centre.position .= [3.0, -0.4, 4.2]
        state.centre.velocity .= [0.03, -0.09, 0.01]
        LilGuys.Centres.calc_centre!(state, snap)
        cen = state.centre

        @test cen.position ≈ [0.,0.,0.] atol=0.5
        @test cen.velocity ≈ [0.,0.,0.] atol=0.005


        snap1 = lguys.deepcopy(snap)
        delta_x = [82.3, -100.5, 391.]
        delta_v = [0.292, -0.051, 0.130]
        snap1.positions .+= delta_x
        snap1.velocities .+= delta_v

        state = lguys.SS_State(snap1)
        LilGuys.Centres.calc_centre!(state, snap1)
        cen = state.centre
        @test cen.position ≈ delta_x atol=0.5
        @test cen.velocity ≈ delta_v atol=0.005
        @test cen.position_err ≈ dx rtol=0.5
        @test cen.velocity_err ≈ dv rtol=0.5
    end

    @testset "stopping criteria rmax" begin
        state = lguys.SS_State(snap, r_max=3)
        @test_logs (:info, r".*completed with status r") match_mode=:any LilGuys.Centres.calc_centre!(state, snap)
        r_max = maximum(radii(snap.positions[:, state.filt], state.centre.position))
        @test r_max ≈ 3 rtol=0.05

        state = lguys.SS_State(snap, r_max=8.5)
        @test_logs (:info, r".*completed with status r") match_mode=:any LilGuys.Centres.calc_centre!(state, snap)
        r_max = maximum(radii(snap.positions[:, state.filt], state.centre.position))
        @test r_max ≈ 8.5 rtol=0.05
    end

    @testset "stopping criteria Nmin" begin
        state = lguys.SS_State(snap, N_min = 300, dx_atol=0, r_factor=0.95)
        @test_logs (:info, r".*completed with status N") match_mode=:any  LilGuys.Centres.calc_centre!(state, snap)
        @test sum(state.filt) ≈ 300 rtol=0.1

        state = lguys.SS_State(snap, N_min = 765, dx_atol=0, r_factor=0.95)
        @test_logs (:info, r".*completed with status N") match_mode=:any LilGuys.Centres.calc_centre!(state, snap) 
        @test sum(state.filt) ≈ 765 rtol=0.1
    end

    @testset "stopping criteria itermax" begin
        state = lguys.SS_State(snap, dx_atol=1e-12, r_factor=1, itermax=10)
        test_logger = TestLogger()

        with_logger(test_logger) do
            LilGuys.Centres.calc_centre!(state, snap)
        end
        @test occursin("itermax", test_logger.logs[end-1].message)
        @test occursin("i=10", test_logger.logs[end-2].message)

        state = lguys.SS_State(snap, dx_atol=0, r_factor=1, itermax=27)
        test_logger = TestLogger()

        with_logger(test_logger) do
            LilGuys.Centres.calc_centre!(state, snap)
        end
        @test occursin("itermax", test_logger.logs[end-1].message)
        @test occursin("i=27", test_logger.logs[end-2].message)
    end


    @testset "stopping criteria dx_rel" begin
        state = lguys.SS_State(snap, dx_atol=1e-12, dx_rtol=0.1)
        @test_logs (:info, r".*completed with status dx_rel") match_mode=:any LilGuys.Centres.calc_centre!(state, snap)
    end

    @testset "stopping criteria dx" begin
        state = lguys.SS_State(snap, dx_atol=0.1, dx_rtol=0.0001)
        @test_logs (:info, r".*completed with status dx") match_mode=:any LilGuys.Centres.calc_centre!(state, snap)
    end

    @testset "stopping criteria dN" begin
        state = lguys.SS_State(snap, dx_atol=0.00001, dN_min=10, r_factor=0.92)
        @test_logs (:info, r".*completed with status dN") (:info, r".*, dN=[0-9],.*")  match_mode=:any LilGuys.Centres.calc_centre!(state, snap)
    end
end




@testset "multimodal" begin
    halo1 = LilGuys.TruncNFW(r_s=2.0,  M_s=0.6, trunc=10)
    halo2 = LilGuys.TruncNFW(r_s=1.0,  M_s=0.1, trunc=10)
    snap1 = LilGuys.sample_potential(halo1, 10000)
    snap2 = LilGuys.sample_potential(halo2, 1000)
    x1 = [0.5, 0.8, 0.12]
    x2 = [-0.2, 1.5, -1.2]
    v1 = [0.05, -0.08, 0.18]
    v2 = [0.02, 0.01, -0.03]

    snap1.positions .+= x1
    snap1.velocities .+= v1
    snap2.positions .+= x2
    snap2.velocities .+= v2

    snap_combined = LilGuys.Snapshot(hcat(snap1.positions, snap2.positions), hcat(snap1.velocities, snap2.velocities), vcat(snap1.masses, snap2.masses))
    snap_combined.potential = LilGuys.potential_spherical_discrete(snap_combined)

    state = lguys.SS_State(snap_combined, x0=x2)
    LilGuys.Centres.calc_centre!(state, snap_combined)
    cen = state.centre

    @test cen.position ≈ x1 atol=0.4
    @test cen.velocity ≈ v1 atol=0.05 # higher because snapshots overlap
end


@testset "pathological" begin
    # no particles, nans, infs
    @test false broken = true 
end


@testset "tidal" begin
    nfw = LilGuys.TruncNFW(r_s=5.5,  M_s=0.8, trunc=12)
    snap = LilGuys.sample_potential(nfw, 8259)

end


@testset "other codes" begin
    # make notebooks to compare against other codes
    @test false broken = true 
end

@testset "runtime scaling" begin
    # make notebooks to time e.g. 1e4-1e6 particles
    @test false broken = true 
end

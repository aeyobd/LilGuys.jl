
using LinearAlgebra: ×
import DensityEstimators: bins_min_width_equal_number



@testset "v_rad" begin
    @testset "basic" begin

        # simple directions
        pos = [ [-1,0,0]     [0, 2, 0]     [1/2, 3/2, -1] [π, √π, π^π]]
        vel = [ [π, 2, -2]  [0, 2.1, 0]   [0, 0, 0]      [0,0,0]]

        @test lguys.radial_velocities(pos, vel) ≈ [-π, 2.1, 0, 0]

        # easy angles
        pos = [ [3, 4, 0] [0, 12/3, -5/3]]
        vel = [ [2, 0, 0] [1, 0, 0.1]]

        @test lguys.radial_velocities(pos, vel) ≈ [2*3/5, 0.1 * -5/13]
    end


    @testset "random" begin
        N = 100
        pos = randn(3, N)

        # random orthoganal
        y = randn(3, N)
        vel = rand(N)' .* hcat([y[:, i] × pos[:, i] for i in 1:N]...)

        @test lguys.radial_velocities(pos, vel) ≈ zeros(N) atol=1e-13

        # random parallel
        v = 10 .^ randn(N)
        vel = v' .* pos ./ lguys.radii(pos)'

        @test lguys.radial_velocities(pos, vel) ≈ v

        pos = pos .* (10 .^ randn(N))'
        @test lguys.radial_velocities(pos, vel) ≈ v

        # antiparallel
        vel = -vel
        @test lguys.radial_velocities(pos, vel) ≈ -v
    end
end



@testset "v_circ_max_model" begin
    param = (1, 1)
    halo = lguys.NFW(v_circ_max=1.0, r_circ_max=1.0)

    @test lguys._v_circ_max_model(1, param) ≈ 1.0
    @test lguys._v_circ_max_model(1.01, param) < 1
    @test lguys._v_circ_max_model(0.99, param) < 1
    @test lguys._v_circ_max_model(-1, param) === NaN

    param = (2.23, √2)
    halo = lguys.NFW(v_circ_max=param[2], r_circ_max=param[1])
    @test lguys._v_circ_max_model(param[1], param) ≈ param[2]

    x = LinRange(0.1, 10, 100)
    y = lguys._v_circ_max_model(x, param)
    y2 = lguys.v_circ.(halo, x)

    @test y ≈ y2

    @test lguys._v_circ_max_model(2.23, (NaN, NaN)) === NaN
    @test lguys._v_circ_max_model(2.23, (-1, 2.4)) === NaN
    @test lguys._v_circ_max_model(π, (π, -2.4)) ≈ -2.4
end


@testset "fit_v_r_circ_max" begin
    @testset "simple cases" begin
        halo = lguys.NFW(v_circ_max=1.0, r_circ_max=1.0)
        r = 10 .^ LinRange(-1, 1, 100)
        v = lguys.v_circ.(halo, r)

        fit = lguys.fit_v_r_circ_max(r, v)
        @test fit.r_circ_max ≈ 1.0 
        @test fit.v_circ_max ≈ 1.0 
        @test fit.converged


        halo = lguys.NFW(v_circ_max=π/2, r_circ_max=√π)
        r = 10 .^ LinRange(-1, 1, 100)
        v = lguys.v_circ.(halo, r)

        fit = lguys.fit_v_r_circ_max(r, v)
        @test fit.r_circ_max ≈ √π
        @test fit.v_circ_max ≈ π/2
        @test fit.converged
    end

    @testset "cannot converge" begin
        r = 10 .^ LinRange(-1, 1, 100)
        v = fill(2.23, 100)

        fit = lguys.fit_v_r_circ_max(r, v)
        @test !fit.converged
        @test fit.v_circ_max ≈ 2.23
    end

    @testset "no data" begin
        fit = lguys.fit_v_r_circ_max([1.], [1.])
        @test !fit.converged
        @test fit.v_circ_max === NaN
        @test fit.r_circ_max === NaN
    end

    @testset "exceptions" begin
        r = [1.0, 2.0, 3.0]
        v = [1.0, 2.0]

        @test_throws DimensionMismatch lguys.fit_v_r_circ_max(r, v)

    end
end



@testset "v_circ" begin
end


@testset "ρ_hist" begin

end


@testset "MassProfile (integration)" begin
    N = 30_000
    M_s = 2
    r_s = 5

    halo = lguys.TruncNFW(M_s=M_s, r_s=r_s, trunc=100)
    M_0 = lguys.mass(halo)

    ρ(r) = lguys.density(halo, r)

    r = lguys.sample_density(ρ, N, log_r=LinRange(-5, 5, 10_000))

    mass = M_0/N  * (1 .+ 0.0randn(N))
    M = sum(mass)

    snap = lguys.Snapshot(positions=r' .* lguys.rand_unit(N), velocities=zeros(3, N), masses=mass, index=1:N, header=lguys.make_default_header(1, N), potential=-ones(N))

    radii = lguys.radii(snap)
    bins = lguys.Interface.bins_both(log10.(radii), nothing, bin_width=0.05, num_per_bin=100)

    profile = lguys.MassProfile3D(snap, bins=bins)

    @test profile.N_bound ≈ N
    @test profile.M_in[end] ≈ M

    @test lguys.mass(halo) ≈ M rtol=1e-2

    r = 10 .^ profile.log_r[2:end-1]
    ρ_exp = lguys.density.(halo, r)
    @test_χ2 profile.rho[2:end-1] profile.rho_err[2:end-1] ρ_exp

    r = 10 .^ profile.log_r_bins[2:end]
    M_exp = lguys.mass.(halo, r)
    @test_χ2 profile.M_in profile.M_in_err M_exp

    # errors tend to be overestimated here...
    v_circ_exp = lguys.v_circ.(halo, profile.r_circ)
    @test_χ2 profile.v_circ profile.v_circ_err v_circ_exp

end



@testset "to_sky" begin
    @testset "inverse" begin
        N = 100
        snap = lguys.Snapshot(positions=100randn(3, N), velocities=1randn(3, N), masses=ones(N), index=1:N, header=lguys.make_default_header(1, N), potential=-ones(N))

        sky = lguys.to_sky(snap)
        gc = lguys.transform.(lguys.Galactocentric, sky)

        positions = [[g.x, g.y, g.z] for g in gc]
        velocities = [[g.v_x, g.v_y, g.v_z] for g in gc]
        positions = hcat(positions...)
        velocities = hcat(velocities...)

        velocities ./= V2KMS

        @test positions ≈ snap.positions
        @test velocities ≈ snap.velocities
    end
end


@testset "to_gaia" begin
    N = 100
    pos = 100randn(3, N)
    vel = 1randn(3, N)
    masses = ones(N)
    weights = rand(N)

    snap = lguys.Snapshot(positions=pos, velocities=vel, masses=masses, weights=weights,
        index=1:N, header=lguys.make_default_header(1, N), potential=-ones(N))

    gaia = lguys.to_gaia(snap, add_centre=false)

    @testset "integration" begin
        obs = lguys.to_sky(snap) 

        @test size(gaia, 1) == N 
        @test gaia.weights ≈ weights
        @test ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"] ⊆ names(gaia)

        @test gaia.distance ≈ [o.distance for o in obs]
        @test gaia.ra ≈ [o.ra for o in obs]

        # test chi2, r_ell, etc
        #
        #
    end

    @testset "add centre" begin
        # test add centre
        gaia_cen = lguys.to_gaia(snap, add_centre=true)
        @test size(gaia_cen, 1) == N + 1
        @test gaia_cen.weights[1] == 0
        @test gaia_cen.distance[1] == 0 broken=true

    end
    
    @testset "filters" begin
        
    end
end


@testset "to_frame" begin
end




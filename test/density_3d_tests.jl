
using LinearAlgebra: ×
import DensityEstimators: bins_min_width_equal_number


@testset "calc_ρ_from_hist" begin
    @testset "simple cases" begin
        masses = [1, 0, π]
        bins = [0, 1, 2, 5]

        ρ = lguys.calc_ρ_from_hist(bins, masses)
        ρ_exp = [3/4π, 0, 3/4 * 1/117]

        @test ρ ≈ ρ_exp 
    end

    @testset "errors" begin
        masses = Int[]
        bins = [1, 2]

        @test_throws DimensionMismatch lguys.calc_ρ_from_hist(masses, bins)

        masses = [1.0, 2.0]
        bins = [1, 2]

        @test_throws DimensionMismatch lguys.calc_ρ_from_hist(masses, bins)
    end

    @testset "mass conservation" begin
        N = 1000
        x = randn(N)
        bins, masses = lguys.histogram(x)

        ρ = lguys.calc_ρ_from_hist(bins, masses)
        M_per_shell = 4π/3*diff(bins .^ 3) .* ρ
        M_tot = sum(masses)

        @test M_tot ≈ N rtol=1e-5
    end
end



@testset "calc_v_rad" begin
    @testset "basic" begin

        # simple directions
        pos = [ [-1,0,0]     [0, 2, 0]     [1/2, 3/2, -1] [π, √π, π^π]]
        vel = [ [π, 2, -2]  [0, 2.1, 0]   [0, 0, 0]      [0,0,0]]

        @test lguys.calc_v_rad(pos, vel) ≈ [-π, 2.1, 0, 0]

        # easy angles
        pos = [ [3, 4, 0] [0, 12/3, -5/3]]
        vel = [ [2, 0, 0] [1, 0, 0.1]]

        @test lguys.calc_v_rad(pos, vel) ≈ [2*3/5, 0.1 * -5/13]
    end


    @testset "random" begin
        N = 100
        pos = randn(3, N)

        # random orthoganal
        y = randn(3, N)
        vel = rand(N)' .* hcat([y[:, i] × pos[:, i] for i in 1:N]...)

        @test lguys.calc_v_rad(pos, vel) ≈ zeros(N) atol=1e-13

        # random parallel
        v = 10 .^ randn(N)
        vel = v' .* pos ./ lguys.calc_r(pos)'

        @test lguys.calc_v_rad(pos, vel) ≈ v

        pos = pos .* (10 .^ randn(N))'
        @test lguys.calc_v_rad(pos, vel) ≈ v

        # antiparallel
        vel = -vel
        @test lguys.calc_v_rad(pos, vel) ≈ -v
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
    y2 = lguys.calc_v_circ.(halo, x)

    @test y ≈ y2

    @test lguys._v_circ_max_model(2.23, (NaN, NaN)) === NaN
    @test lguys._v_circ_max_model(2.23, (-1, 2.4)) === NaN
    @test lguys._v_circ_max_model(π, (π, -2.4)) ≈ -2.4
end


@testset "fit_v_r_circ_max" begin
    @testset "simple cases" begin
        halo = lguys.NFW(v_circ_max=1.0, r_circ_max=1.0)
        r = 10 .^ LinRange(-1, 1, 100)
        v = lguys.calc_v_circ.(halo, r)

        fit = lguys.fit_v_r_circ_max(r, v)
        @test fit.r_circ_max ≈ 1.0
        @test fit.v_circ_max ≈ 1.0
        @test fit.converged


        halo = lguys.NFW(v_circ_max=π/2, r_circ_max=√π)
        r = 10 .^ LinRange(-1, 1, 100)
        v = lguys.calc_v_circ.(halo, r)

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

    @testset "exceptions" begin
        r = [1.0, 2.0, 3.0]
        v = [1.0, 2.0]

        @test_throws DimensionMismatch lguys.fit_v_r_circ_max(r, v)

        @test_throws ArgumentError lguys.fit_v_r_circ_max([1.], [1.])
    end
end



@testset "calc_v_circ" begin
end


@testset "calc_ρ_hist" begin

end

@testset "calc_M_in" begin
    @testset "simple_cases" begin
        N = 10
        snap = lguys.Snapshot(positions=zeros(3, N), velocities=zeros(3, N), masses=ones(N), index=1:N, header=lguys.make_default_header(1, N), Φs=-ones(N))

        @test lguys.calc_M_in(snap, 0.01) == N


        r = collect(1:10) .- 0.01
        snap = lguys.Snapshot(positions=r' .* lguys.rand_unit(N), velocities=zeros(3, N), masses=ones(N), index=1:N, header=lguys.make_default_header(1, N), Φs=-ones(N))

        expected = collect(1:10)
        actual = [lguys.calc_M_in(snap, i) for i in 1:N]
        @test expected ≈ actual
    end
end



@testset "MassProfile (integration)" begin
    N = 30_000
    M_s = 2
    r_s = 5

    halo = lguys.TruncNFW(M_s=M_s, r_s=r_s, trunc=100)
    M_0 = lguys.calc_M_tot(halo)

    ρ(r) = lguys.calc_ρ(halo, r)

    r = lguys.sample_ρ(ρ, N, log_r=LinRange(-5, 5, 10_000))

    mass = M_0/N  * (1 .+ 0.0randn(N))
    M = sum(mass)

    snap = lguys.Snapshot(positions=r' .* lguys.rand_unit(N), velocities=zeros(3, N), masses=mass, index=1:N, header=lguys.make_default_header(1, N), Φs=-ones(N))

    radii = lguys.calc_r(snap)
    bins = bins_min_width_equal_number(log10.(radii), dx_min=0.05, N_per_bin_min=100)#[2:end-1]

    profile = lguys.MassProfile3D(snap, bins=bins)

    @test profile.N_bound ≈ N
    @test profile.M_in[end] ≈ M

    @test lguys.get_M_tot(halo) ≈ M rtol=1e-2

    r = 10 .^ profile.log_r[2:end-1]
    ρ_exp = lguys.calc_ρ.(halo, r)
    @test_χ2 profile.rho[2:end-1] profile.rho_err[2:end-1] ρ_exp

    r = 10 .^ profile.log_r_bins[2:end]
    M_exp = lguys.calc_M.(halo, r)
    @test_χ2 profile.M_in profile.M_in_err M_exp

    # errors tend to be overestimated here...
    v_circ_exp = lguys.calc_v_circ.(halo, profile.r_circ)
    @test_χ2 profile.v_circ profile.v_circ_err v_circ_exp

end



@testset "to_sky" begin
    @testset "inverse" begin
        N = 100
        snap = lguys.Snapshot(positions=100randn(3, N), velocities=1randn(3, N), masses=ones(N), index=1:N, header=lguys.make_default_header(1, N), Φs=-ones(N))

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
        index=1:N, header=lguys.make_default_header(1, N), Φs=-ones(N))

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




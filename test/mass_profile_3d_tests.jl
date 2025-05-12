
using LinearAlgebra: ×
import DensityEstimators: bins_min_width_equal_number


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

    profile = lguys.MassProfile(snap, bins=bins)
    sc = lguys.MassScalars(snap, profile)

    @test sc.N_bound ≈ N
    @test profile.M_in[end].middle ≈ M rtol=1e-2

    @test lguys.mass(halo) ≈ M rtol=1e-2

    r = lguys.radii(profile)
    M_exp = lguys.mass.(halo, r)
    @test_χ2 profile.M_in M_exp

    # errors tend to be overestimated here...
    v_circ_exp = lguys.v_circ.(halo, lguys.radii(profile))
    @test_χ2 lguys.circular_velocity(profile) v_circ_exp

end





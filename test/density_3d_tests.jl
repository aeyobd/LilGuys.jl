
using LinearAlgebra: ×
import DensityEstimators: bins_min_width_equal_number


@testset "density_from_hist" begin
    @testset "simple cases" begin
        masses = [1, 0, π]
        bins = [0, 1, 2, 5]

        ρ = lguys.density_from_hist(bins, masses)
        ρ_exp = [3/4π, 0, 3/4 * 1/117]

        @test ρ ≈ ρ_exp 
    end

    @testset "errors" begin
        masses = Int[]
        bins = [1, 2]

        @test_throws DimensionMismatch lguys.density_from_hist(masses, bins)

        masses = [1.0, 2.0]
        bins = [1, 2]

        @test_throws DimensionMismatch lguys.density_from_hist(masses, bins)
    end

    @testset "mass conservation" begin
        N = 1000
        x = randn(N)
        bins, masses = lguys.histogram(x)

        ρ = lguys.density_from_hist(bins, masses)
        M_per_shell = 4π/3*diff(bins .^ 3) .* ρ
        M_tot = sum(masses)

        @test M_tot ≈ N rtol=1e-5
    end
end



@testset "DensityProfile3D (integration)" begin
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

    profile = lguys.DensityProfile3D(snap, bins=bins)

    @test sum(profile.counts) ≈ N

    r = 10 .^ profile.log_r[2:end-1]
    ρ_exp = lguys.density.(halo, r)
    @test_χ2 profile.rho[2:end-1] profile.rho_err[2:end-1] ρ_exp
end



@testset "radial_velocities" begin
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


@testset "spherical_unit_vectors" begin
    @test false broken = true
end


@testset "to_spherical_velocities" begin
    @test false broken = true
end


@testset "β_profile" begin
    @test false broken = true
end

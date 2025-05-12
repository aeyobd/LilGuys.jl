@testset "σv" begin
    N = 10000
    σ_exp = 2.1
    positions = randn(3, N)
    velocities = σ_exp * randn(3, N)
    masses = rand(N)
    weights = rand(N)

    snap = lguys.Snapshot(positions, velocities, masses, weights=weights)

    @test lguys.σv_1d(snap) ≈ σ_exp rtol=1e-2
    @test lguys.σv_x(snap) ≈ σ_exp rtol=3e-2


    r_max = 0.5
    filt = lguys.radii(positions) .< r_max

    @test sum(.!filt) > 0
    velocities[:, .!filt] .*= 3 

    @test lguys.σv_1d(snap, r_max=r_max) ≈ σ_exp rtol=1e-2
    @test lguys.σv_1d(snap, r_max=2) > σ_exp * 1.5

    @test lguys.σv_x(snap, r_max=r_max) ≈ σ_exp rtol=3e-2
    @test lguys.σv_x(snap, r_max=2) > σ_exp * 1.5
end



@testset "break_radius"  begin
    @test lguys.break_radius(1.0, 1.0) ≈ 0.55
    @test lguys.break_radius(0.5, 1.0) ≈ 0.55 * 0.5
    @test lguys.break_radius(1.2, 0.5) ≈ 0.55 * 0.6
end



@testset "stellar profile 3D (integration)" begin
    ρ(r) = exp(-r)/8π
    N = 10_000
    r = LilGuys.sample_density(ρ, N, log_r=LinRange(-5, 5, 1000))

    positions = r' .* lguys.rand_unit(N)
    σv_exp = 0.055

    velocities = σv_exp .* randn(3, N)
    masses = rand(N)

    weights = 1.75 .+ 0.5 * rand(N) # mean of 2.0
    weights ./= N

    Mtot = sum(weights)

    snap = lguys.Snapshot(positions, velocities, masses, weights=weights)

    bins = 20
    prof = lguys.DensityProfile(snap, snap.weights, bins=bins)
    @test prof isa lguys.DensityProfile

    @testset "properties" begin
        ρ_exp = @. Mtot .* ρ(10 ^ prof.log_r)

        @test_χ2 prof.rho ρ_exp

        @test length(prof.log_r) == bins 
        @test length(prof.log_r_bins) == bins  + 1
    end


    @testset "scaling" begin
        m_scale = 1.5
        r_scale = 0.95
        m_scale_pot = 1.1

        v_scale = sqrt(m_scale_pot / r_scale)
        snap_scaled = lguys.Snapshot(positions * r_scale, velocities * v_scale, masses, weights=weights * m_scale)
        prof_scaled = lguys.DensityProfile(snap_scaled, snap_scaled.weights, bins=bins)
        prof_rescaled = lguys.scale(prof, r_scale, m_scale, m_scale_pot)

        for k in propertynames(prof_rescaled)
            if k != :annotations
                if eltype(getproperty(prof_rescaled, k)) <: Measurement
                    @test middle.(getproperty(prof_rescaled, k)) ≈ middle.(getproperty(prof_scaled, k)) nans=true
                else
                    @test getproperty(prof_rescaled, k) ≈ getproperty(prof_scaled, k) rtol=1e-2 nans=true
                end
            end
        end
    end



    @testset "stellar scalars" begin 
        s = lguys.StellarScalars(snap, 
            r_max = 3, delta_t=2.0
           )

        @test lguys.time_since_peri(s) == 2.0
        @test lguys.break_radius(s) ≈ 2.0 * 0.55 * σv_exp atol=1e-2
        @test lguys.velocity_dispersion(s) ≈ σv_exp rtol=3e-2
        @test s.r_max_sigma ≈ 3.0
        @test s.time ≈ snap.time
        @test s.bound_mass ≈ 2.0 rtol=1e-3
    end
end



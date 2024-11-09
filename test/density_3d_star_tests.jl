@testset "initialization" begin

end


@testset "calc_σv" begin
    N = 10000
    σ_exp = 2.1
    positions = randn(3, N)
    velocities = σ_exp * randn(3, N)
    masses = rand(N)
    weights = rand(N)

    snap = lguys.Snapshot(positions, velocities, masses, weights=weights)

    @test lguys.calc_σv_1d(snap) ≈ σ_exp rtol=1e-2
    @test lguys.calc_σv_x(snap) ≈ σ_exp rtol=3e-2


    r_max = 0.5
    filt = lguys.calc_r(positions) .< r_max

    @test sum(.!filt) > 0
    velocities[:, .!filt] .*= 3 

    @test lguys.calc_σv_1d(snap, r_max=r_max) ≈ σ_exp rtol=1e-2
    @test lguys.calc_σv_1d(snap, r_max=2) > σ_exp * 1.5

    @test lguys.calc_σv_x(snap, r_max=r_max) ≈ σ_exp rtol=3e-2
    @test lguys.calc_σv_x(snap, r_max=2) > σ_exp * 1.5
end



@testset "calc_break_radius"  begin
    @test lguys.calc_break_radius(1.0, 1.0) ≈ 0.55
    @test lguys.calc_break_radius(0.5, 1.0) ≈ 0.55 * 0.5
    @test lguys.calc_break_radius(1.2, 0.5) ≈ 0.55 * 0.6
end



@testset "stellar profile 3D (integration)" begin
    ρ(r) = exp(-r)/8π
    N = 10_000
    r = LilGuys.sample_ρ(ρ, N, log_r=LinRange(-5, 5, 1000))

    positions = r' .* lguys.rand_unit(N)
    σv_exp = 0.055

    velocities = σv_exp .* randn(3, N)
    masses = rand(N)

    weights = 1.75 .+ 0.5 * rand(N) # mean of 2.0
    weights ./= N

    Mtot = sum(weights)

    snap = lguys.Snapshot(positions, velocities, masses, weights=weights)

    bins = 20
    prof = lguys.StellarProfile3D(snap, bins=bins)
    @test prof isa lguys.StellarProfile3D

    @testset "properties" begin
        @test sum(prof.mass_in_shell) ≈ Mtot rtol=1e-2
        ρ_exp = @. Mtot .* ρ(10 ^ prof.log_r)

        @test_χ2 prof.rho prof.rho_err ρ_exp

        @test prof.sigma_vx ≈ σv_exp rtol=1e-2

        @test prof.quantiles[3:5] ≈ [0.1, 0.5, 0.9]
        @test prof.r_quantile[3:5] ≈ [1.102, 2.6740, 5.3223] rtol=1e-2
        # quantiles test
        # integral of the profile is M * 1/2 * (2 - (x^2 + 2x + 2) * exp(-x))
        @test length(prof.log_r) == bins 
        @test length(prof.log_r_bins) == bins  + 1
    end


    @testset "scaling" begin
        m_scale = 1.5
        r_scale = 0.95
        v_scale = 1.1

        snap_scaled = lguys.Snapshot(positions * r_scale, velocities * v_scale, masses, weights=weights * m_scale)
        prof_scaled = lguys.StellarProfile3D(snap_scaled, bins=bins)
        prof_rescaled = lguys.scale(prof, r_scale, v_scale, m_scale)

        for k in propertynames(prof_rescaled)
            @test getproperty(prof_rescaled, k) ≈ getproperty(prof_scaled, k) rtol=1e-2 nans=true
        end
    end

    @testset "arguments" begin

        prof = lguys.StellarProfile3D(snap, bins=bins, 
            r_max = 3, quantiles = [0.3, 0.5, 0.7], delta_t=2.0
           )

        @test prof.quantiles == [0.3, 0.5, 0.7]
        @test prof.delta_t == 2.0
        @test prof.r_break ≈ 2.0 * 0.55 * σv_exp atol=1e-2
        @test prof.sigma_vx ≈ lguys.calc_σv_1d(snap, r_max=3) rtol=1e-8
    end


end

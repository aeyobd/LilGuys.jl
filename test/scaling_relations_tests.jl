@testset "fattahi2018" begin
    @test log10( LilGuys.M_s_from_vel_fattahi(100 / V2KMS) ) ≈ -0.5 atol=0.1
    @test log10( LilGuys.M_s_from_vel_fattahi(50 / V2KMS) ) ≈ -2 atol=0.1
    @test log10( LilGuys.M_s_from_vel_fattahi(30 / V2KMS) ) ≈ -3.8 atol=0.1
    @test log10( LilGuys.M_s_from_vel_fattahi(18 / V2KMS) ) ≈ -8 atol=0.1
end


@testset "fattahi2018 inverse" begin
    log_Ms = LinRange(-10, 0, 100)
    vels = LilGuys.vel_from_M_s_fattahi.(10 .^ log_Ms)
    log_Ms_recovered = log10.( LilGuys.M_s_from_vel_fattahi.(vels) )

    @test log_Ms_recovered ≈ log_Ms rtol=0.001
end


@testset "EN21" begin
    @test LilGuys.v_circ_EN21(1) ≈ 1.0
    @test LilGuys.v_circ_EN21(0.1) ≈ 10 ^ -0.54 rtol=0.03
    @test LilGuys.v_circ_EN21(10^-0.4) ≈ 10 ^ -0.16 rtol=0.03
end

@testset "EN21_tidal_track" begin
    r0 = 2.32
    v0 = 0.153
    rcm, vcm = LilGuys.EN21_tidal_track(r0, v0)
    @test rcm[1] ≈ r0 rtol=1e-2
    @test vcm[1] ≈ v0 rtol=1e-2

    @test vcm[end] / v0 ≈ LilGuys.v_circ_EN21(rcm[end]/r0) rtol=1e-7
end


@testset "Ludlow" begin
    solve_rmax = LilGuys.Ludlow.solve_rmax
    c_ludlow = LilGuys.Ludlow.c_ludlow

    # tests from Asya's papers
    # Crater II (Borukhovetskaya et al. 2022)
    @test solve_rmax(25.9/V2KMS) ≈ 4.7 atol = 0.05

    # Fornax (Borukhovetskaya et al. 2022)
    @test solve_rmax(39.6/V2KMS) ≈ 8.0 atol = 0.05
end



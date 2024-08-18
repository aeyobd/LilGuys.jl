
@testset "A nfw" begin
    @test lguys.A_NFW(1) ≈ 0.1931471805599453
    @test lguys.A_NFW(5.6) ≈ 1 rtol=1e-1
    @test lguys.A_NFW(Inf) === Inf
    @test lguys.A_NFW(0.) == 0 

    @test_throws DomainError lguys.A_NFW(-2) 
    @test_throws DomainError lguys.A_NFW(-Inf) 
    @test_throws DomainError lguys.A_NFW(-0.000000001)
end


@testset "NFW creation" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test nfw.r_s == 1
    @test nfw.M_s == 1
    @test nfw.c > 0
end


@testset "NFW alternate creations" begin
    halo = lguys.NFW(r_s=8.5, M_s=1.1)
    @test halo.r_s == 8.5
    @test halo.M_s == 1.1

    M200 = lguys.calc_M200(halo)
    R200 = lguys.calc_R200(halo)
    v_circ_max = lguys.calc_v_circ_max(halo)
    r_circ_max = lguys.calc_r_circ_max(halo)

    halo2 = lguys.NFW(M200=M200, r_s=halo.r_s)
    halo3 = lguys.NFW(v_circ_max=v_circ_max, r_circ_max=r_circ_max)
    halo4 = lguys.NFW(M200=M200, c=halo.c)

    for alt_halo in [halo2, halo3, halo4]
        @test alt_halo.r_s ≈ halo.r_s
        @test alt_halo.M_s ≈ halo.M_s
        @test alt_halo.c ≈ halo.c
    end
end

@testset "NFW creation error" begin
    @test_throws ArgumentError lguys.NFW(r_s=1, M_s=1, M200=1)
    @test_throws ArgumentError lguys.NFW(r_s=1, M_s=1, R200=1)
    @test_throws ArgumentError lguys.NFW()
    @test_throws ArgumentError lguys.NFW(r_s=1)
    @test_throws ArgumentError lguys.NFW(M200=Nothing, R200=nothing)
end


@testset "NFW density" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test lguys.calc_ρ(nfw, 0) === Inf
    @test lguys.calc_ρ(nfw, 1) ≈ 1/(16π)
end



@testset "NFW mass" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test lguys.calc_M(nfw, 0) == 0
    @test lguys.calc_M(nfw, Inf) === Inf
end



@testset "NFW potential" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test lguys.calc_Φ(nfw, 0) == -1
    @test lguys.calc_Φ(nfw, 1e-10) ≈ -1
    @test lguys.calc_Φ(nfw, 1) ≈ -log(2)
    @test lguys.calc_Φ(nfw, 1e10) ≈ 0 atol=1e-8
    @test lguys.calc_Φ(nfw, Inf) == 0.
    @test lguys.calc_Φ(nfw, NaN) === NaN

    @test_throws DomainError lguys.calc_Φ(nfw, -1)
end


@testset "v circ max" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    r_max = lguys.calc_r_circ_max(nfw)
    v_max = lguys.calc_v_circ_max(nfw)

    # is local maximum and agrees roughly
    @test lguys.calc_v_circ(nfw, r_max) ≈ v_max
    @test lguys.calc_v_circ(nfw, r_max * (1.01)) < v_max
    @test lguys.calc_v_circ(nfw, r_max * (0.99)) < v_max
end


@testset "NFW scaling" begin
    nfw = lguys.NFW(r_s=1, M_s=1)

    radii = 10 .^ range(-2, 2, length=100)
    M_exp = lguys.calc_M.(nfw, radii)
    ρ_exp = lguys.calc_ρ.(nfw, radii)
    Φ_exp = lguys.calc_Φ.(nfw, radii)

    for _ in 1:100
        M_s = 10 ^ (randn()/4)
        r_s = 10 ^ (randn()/4 + 1)

        halo = lguys.NFW(r_s=r_s, M_s=M_s)

        r = radii * r_s
        M = lguys.calc_M.(halo, r)
        ρ = lguys.calc_ρ.(halo, r)
        Φ = lguys.calc_Φ.(halo, r)

        @test M ≈ M_exp * M_s
        @test ρ ≈ ρ_exp * M_s / r_s^3
        @test Φ ≈ Φ_exp * M_s / r_s

    end
end


@testset "NFW 200" begin
    halo = lguys.NFW(r_s=10.53, M_s=π*3/2)
    M200 = lguys.calc_M200(halo)
    R200 = lguys.calc_R200(halo)

    ρ_mean = M200 / (4/3 * π * R200^3)
    @test ρ_mean ≈ 200*lguys.ρ_crit
    @test lguys.calc_M(halo, R200) ≈ M200
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


@testset "literature datapoints" begin

    # Fornax (Borukhovetskaya et al. 2022)
    halo = lguys.NFW(v_circ_max = 39.6/V2KMS, r_circ_max = 8.0 / R2KPC)
    @test lguys.calc_M200(halo) ≈ 1.04 atol=0.01
    @test halo.c ≈ 12.5 atol=0.05

    # Crater II (Borukhovetskaya et al. 2022)
    halo = lguys.NFW(v_circ_max = 25.9/V2KMS, r_circ_max = 4.7 / R2KPC)
    @test lguys.calc_M200(halo) ≈ 0.272 atol=0.01
    @test halo.c ≈ 13.6 atol=0.05
end

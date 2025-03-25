
@testset "A nfw" begin
    @test lguys.A_NFW(1) ≈ 0.1931471805599453
    @test lguys.A_NFW(5.6) ≈ 1 rtol=1e-1
    @test lguys.A_NFW(Inf) === Inf
    @test lguys.A_NFW(0.) == 0 

    @test lguys.A_NFW(-2)  === NaN
    @test lguys.A_NFW(-Inf)  === NaN
    @test lguys.A_NFW(-1e-10)  === NaN
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

    M200 = lguys.M200(halo)
    R200 = lguys.R200(halo)
    v_circ_max = lguys.v_circ_max(halo)
    r_circ_max = lguys.r_circ_max(halo)

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
    @test lguys.density(nfw, 0) === Inf
    @test lguys.density(nfw, 1) ≈ 1/(16π)
end



@testset "NFW mass" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test lguys.mass(nfw, 0) == 0
    @test lguys.mass(nfw, Inf) === Inf
end



@testset "NFW potential" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test lguys.potential(nfw, 0) == -1
    @test lguys.potential(nfw, 1e-10) ≈ -1
    @test lguys.potential(nfw, 1) ≈ -log(2)
    @test lguys.potential(nfw, 1e10) ≈ 0 atol=1e-8
    @test lguys.potential(nfw, Inf) == 0.
    @test lguys.potential(nfw, NaN) === NaN

    @test_throws DomainError lguys.potential(nfw, -1)
end


@testset "NFW properties" begin
    nfw = lguys.NFW(r_s=3.9, M_s=0.42)

    v1 = lguys.v_circ_max(nfw)
    r = lguys.solve_r_circ_max(nfw)
    v2 = lguys.v_circ(nfw, r)
    @test v1 ≈ v2

    r = lguys.R200(nfw)
    r1 = lguys.solve_R200(nfw)
    @test r ≈ r1
    M = lguys.M200(nfw)
    @test lguys.mass(nfw, r1) ≈ M

    @test lguys.mass(nfw, r1) / (4/3 * π * r1^3) ≈ 200 * lguys.ρ_crit
end

@testset "v circ max" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    r_max = lguys.r_circ_max(nfw)
    v_max = lguys.v_circ_max(nfw)

    # is local maximum and agrees roughly
    @test lguys.v_circ(nfw, r_max) ≈ v_max
    @test lguys.v_circ(nfw, r_max * (1.01)) < v_max
    @test lguys.v_circ(nfw, r_max * (0.99)) < v_max
end


@testset "NFW scaling" begin
    nfw = lguys.NFW(r_s=1, M_s=1)

    radii = 10 .^ range(-2, 2, length=100)
    M_exp = lguys.mass.(nfw, radii)
    ρ_exp = lguys.density.(nfw, radii)
    Φ_exp = lguys.potential.(nfw, radii)

    for _ in 1:100
        M_s = 10 ^ (randn()/4)
        r_s = 10 ^ (randn()/4 + 1)

        halo = lguys.NFW(r_s=r_s, M_s=M_s)

        r = radii * r_s
        M = lguys.mass.(halo, r)
        ρ = lguys.density.(halo, r)
        Φ = lguys.potential.(halo, r)

        @test M ≈ M_exp * M_s
        @test ρ ≈ ρ_exp * M_s / r_s^3
        @test Φ ≈ Φ_exp * M_s / r_s

    end
end


@testset "NFW 200" begin
    halo = lguys.NFW(r_s=10.53, M_s=π*3/2)
    M200 = lguys.M200(halo)
    R200 = lguys.R200(halo)

    ρ_mean = M200 / (4/3 * π * R200^3)
    @test ρ_mean ≈ 200*lguys.ρ_crit
    @test lguys.mass(halo, R200) ≈ M200
end



@testset "TruncNFW" begin
    @testset "simple" begin
        halo = lguys.TruncNFW(r_s=1, M_s=1, r_t=2)

        @test lguys.density(halo, 0) === Inf
        @test lguys.density(halo, 1) ≈ 1/4π * 1/(2^2) * exp(-1/2)
        @test lguys.density(halo, 2) ≈ 1/4π * 1/(2*3^2) * exp(-2/2)
    end


    @testset "Mtot" begin
        halo = lguys.TruncNFW(r_s=2.231, M_s=0.952, r_t=6)

        @test lguys.mass(halo, 0) ≈ 0
        @test lguys.mass(halo, 1000) ≈ lguys.mass(halo)
        @test lguys.mass(halo, 10^4) ≈ lguys.mass(halo)

    end

    @testset "M" begin
        profile = lguys.TruncNFW(r_s=1.21, M_s=√π, trunc=5)
        x = 10 .^ LinRange(-1, 2, 10)

        M1 = lguys.mass.(profile, x)
        M2 = lguys.mass_from_density.(profile, x)

        @test M1 ≈ M2 rtol=1e-5
    end

    @testset "Φ" begin
        profile = lguys.TruncNFW(r_s=1.21, M_s=√π, trunc=5)

        x = 10 .^ LinRange(0, 1, 10)
        Φ1 = lguys.potential.(profile, x)
        Φ2 = lguys.potential_from_density.(profile, x)

        @test Φ1 ≈ Φ2 rtol=1e-5
    end
end


@testset "CoredNFW" begin
    @testset "simple" begin
        halo = lguys.CoredNFW(r_s=1, M_s=1, r_t=2, r_c=0.1)

        @test lguys.density(halo, 0) ≈ 10 / (4π)
        @test lguys.density(halo, 1e-8) ≈ 10 / (4π) rtol=1e-5
        @test lguys.density(halo, 1e-6) ≈ 10 / (4π) rtol=1e-3

        @test lguys.density(halo, 1) ≈ 1/4π * 1/(1.1 * 2^2) * exp(-1/2)
        @test lguys.density(halo, 2) ≈ 1/4π * 1/(2.1 *3^2) * exp(-2/2)

        @test lguys.density(halo, 20) ≈ 0 atol=1e-6
    end

    @testset "M" begin
        profile = lguys.CoredNFW(r_s=1.21, M_s=√π, r_t=5, r_c=0.3)
        x = 10 .^ LinRange(-1, 2, 10)

        M1 = lguys.mass.(profile, x)
        M2 = lguys.mass_from_density.(profile, x)

        @test M1 ≈ M2 rtol=1e-5
    end

    @testset "v circ max" begin
        profile = lguys.CoredNFW(r_s=1.5, M_s=1/π, r_t=6, r_c=0.55)
        r_max = lguys.r_circ_max(profile)
        v_max = lguys.v_circ_max(profile)
        @test lguys.v_circ(profile, r_max) ≈ v_max
        @test lguys.v_circ(profile, r_max * (1.001)) < v_max
        @test lguys.v_circ(profile, r_max * (0.999)) < v_max
    end

    @testset "M200" begin
        profile = lguys.CoredNFW(r_s=2.3, M_s=1/1.5π, r_t=8, r_c=0.03)
        M200 = lguys.M200(profile)
        R200 = lguys.R200(profile)

        ρ_mean = M200 / (4/3 * π * R200^3)
        @test ρ_mean ≈ 200*lguys.ρ_crit
        @test lguys.mass(profile, R200) ≈ M200
    end
end

@testset "literature datapoints" begin

    # Fornax (Borukhovetskaya et al. 2022)
    halo = lguys.NFW(v_circ_max = 39.6/V2KMS, r_circ_max = 8.0 / R2KPC)
    @test lguys.M200(halo) ≈ 1.04 atol=0.01
    @test halo.c ≈ 12.5 atol=0.05

    # Crater II (Borukhovetskaya et al. 2022)
    halo = lguys.NFW(v_circ_max = 25.9/V2KMS, r_circ_max = 4.7 / R2KPC)
    @test lguys.M200(halo) ≈ 0.272 atol=0.01
    @test halo.c ≈ 13.6 atol=0.05
end





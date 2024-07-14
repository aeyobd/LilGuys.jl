
test_profiles = [
    lguys.Exp2D(1, 1),
    lguys.Exp2D(1.3, 0.7124),
    lguys.Exp3D(1.3, 0.7124),
    lguys.Exp3D(1, 1),
    lguys.LogCusp2D(1, 1),
]


@testset "total mass" begin
    for profile in test_profiles
        @test lguys.calc_M(profile, 100.0) ≈ profile.M
    end
end


@testset "total 2D mass" begin
    for profile in test_profiles
        integrand(r) = 2 * π * r * lguys.calc_Σ(profile, r)

        @test lguys.quadgk(integrand, 0, Inf)[1] ≈ profile.M
    end
end


@testset "3d to 2d density" begin

    for profile in test_profiles
        integrand(r, R) = 2 *r* lguys.calc_ρ(profile, r) / sqrt(r^2 - R^2)

        x = 10 .^ LinRange(-2, 0, 100)
        eps = 1e-6
        Σ(R) = lguys.quadgk(r -> integrand(r, R), R*(1+eps), Inf)[1]

        @test lguys.calc_Σ.([profile], x) ≈ Σ.(x) rtol=1e-3
    end
end


@testset "NFW" begin
    prof = lguys.NFW(M_s=1, r_s=1)

    @test lguys.calc_M(prof, 1) ≈ prof.M_s * lguys.A_NFW(1)
end



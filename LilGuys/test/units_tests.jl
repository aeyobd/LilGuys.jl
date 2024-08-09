@testset "Defined values" begin
    @test LilGuys.G == 1
    @test LilGuys.M2MSUN ≈ 1e10
    @test LilGuys.R2KPC ≈ 1

end

@testset "unit consistency" begin
    # validates the last two units: time and velocity
    
    # assuming CODATA 2018 values for G, IAU resolution for MSun, au, year
    # Note: relative error on sqrt(G) is 0.005 which propogates through, so
    # high rtols are okay
    kms_per_kpc_gyr = 0.977792221
    @test V2KMS / kms_per_kpc_gyr ≈ R2KPC / T2GYR rtol=2e-4


    G_physical = 4.300917e-6 # ± 0.000097e-6 kpc^3 Gyr^-2 Msun^-1
    @test G_physical * M2MSUN / R2KPC / V2KMS^2 ≈ 1 rtol=2e-4
end


@testset "arcmin to kpc" begin
    @test LilGuys.arcmin_to_kpc(1, 1) ≈ 0.0002908879 rtol=1e-5
    @test LilGuys.arcmin_to_kpc(40, 86.2) ≈ 1 rtol=3e-3
    @test LilGuys.arcmin_to_kpc(40, 862) ≈ 10 rtol=3e-3
end


@testset "kpc to arcmin" begin
    @test LilGuys.kpc_to_arcmin(1, 1000) ≈ 3.43775 rtol = 1e-5
    @test LilGuys.kpc_to_arcmin(1, 86.2) ≈ 40 rtol = 3e-3
    @test LilGuys.kpc_to_arcmin(0.1, 86.2) ≈ 4 rtol = 3e-3
end


@testset "kpc_to_arcmin inverse" begin
    for (a, b) in [
        (1, 1), (10, 234.2), (1.23, 2344)
       ]
        @test LilGuys.arcmin_to_kpc(LilGuys.kpc_to_arcmin(a, b), b) ≈ a
    end
end

@testset "pm_to_kms" begin
    @test LilGuys.pm_to_kms(1, 1) ≈ 4.74047 rtol=1e-5
    @test LilGuys.pm_to_kms(3, 1) ≈ 3*4.74047 rtol=1e-5
    @test LilGuys.pm_to_kms(3, 0.5) ≈ 3 * 0.5 * 4.74047 rtol=1e-5
end

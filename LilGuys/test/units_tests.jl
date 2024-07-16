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

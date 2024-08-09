
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


@testset "NFW density" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test lguys.calc_ρ(nfw, 0) === Inf
end



@testset "NFW mass" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test lguys.calc_M(nfw, 0) == 0
end



@testset "NFW potential" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    @test lguys.calc_Φ(nfw, 0) === -Inf

end


@testset "v circ max" begin
    nfw = lguys.NFW(r_s=1, M_s=1)
    r_max = lguys.calc_r_circ_max(nfw)
end

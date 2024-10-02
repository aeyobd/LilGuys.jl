
@testset "phase volume" begin
    rs = [0.1, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3]
    vs = zeros(7)
    δr, δv = lguys.phase_volume(rs, vs)


    expected = 0.1
    @test δr ≈ expected rtol=0.1
    @test δv ≈ 0 atol=1e-10


    N = 1000
    k = 10
    δrs = zeros(N)
    δvs = zeros(N)
    for i in 1:N
        rs = sort(sqrt.(rand(k)))
        vs = rand(k)
        δr, δv = lguys.phase_volume(rs, vs)
        δrs[i] = δr
        δvs[i] = δv
    end

    δr = lguys.mean(δrs)
    δv = lguys.mean(δvs)
    expected = N^(-1/3)
    @test δr ≈ expected rtol=0.1
    @test δv ≈ expected rtol=0.1

end


@testset "phase volumes" begin
    N = 6
    snap = uniform_snap(N)

    δr, δv = lguys.phase_volumes(snap, k=4)

    μ_r = lguys.mean(δr)
    σ_r = lguys.std(δr)

    println("r: $μ_r ± $σ_r")
    actual = μ_r * N^(1/3)
    expected = 1.0
    @test actual ≈ expected rtol=0.1

end


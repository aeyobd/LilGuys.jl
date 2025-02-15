import SpecialFunctions: gamma
@testset "phase volume" begin
    rs = [0.1, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3]
    vs = zeros(7)
    k = length(rs)

    δr, δv = lguys.Centres._phase_volume(rs, vs; gamma_ratio=gamma(k+1/3) / gamma(k))


    expected = 0.3 / cbrt(7)
    @test δr ≈ expected rtol=1e-1
    @test δv ≈ 0 atol=1e-10


    N = 30000
    k = 10
    δrs = zeros(N)
    δvs = zeros(N)

    σv = 0.24 

    gamma_ratio = gamma(k+1/3) / gamma(k)
    for i in 1:N
        rs = sort(cbrt.(k * rand(k)))
        vs = σv * abs.(randn(k))
        δr, δv = lguys.Centres._phase_volume(rs, vs, η=1, gamma_ratio=gamma_ratio)
        δrs[i] = δr
        δvs[i] = δv
    end

    δr = lguys.mean(δrs)
    δv = lguys.mean(δvs)
    expected = 1
    @test δr ≈ expected rtol=0.1
    @test δv ≈ σv rtol=0.1

end

function uniform_snap(N; R=1.0, σv=0.36)

    # gaussian distribution
    vs = σv * randn(N)

    rs = cbrt.(R * rand(N))
    positions = rs' .* lguys.rand_unit(N)
    velocities = vs' .* lguys.rand_unit(N)

    return lguys.Snapshot(positions, velocities, 1)
end

@testset "phase volumes" begin
    N = 1000_000
    σv = 0.36
    R = 2.1
    snap = uniform_snap(N; σv=σv, R=R)

    for k in [2, 5, 10, 25]
        δr, δv = lguys.Centres.phase_volumes(snap, k=k)

        filt = lguys.calc_r(snap.positions) .< 0.5*R

        μ_r = lguys.mean(δr[filt])
        σ_r = lguys.std(δr[filt])
        μ_v = lguys.mean(δv)

        actual = μ_r * cbrt(N )
        expected = R  / cbrt(4π/3) # inverse of cbrt of density

        @test actual ≈ expected rtol=0.03

        @test μ_v ≈ σv * √2 rtol=0.1
    end
end


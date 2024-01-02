
@testset "centroid" begin
    x = [1. 5;
         0.5 1.5;
         -1 1;]

    expected = [3., 1, 0]
    cen = lguys.centroid(x)

    @test cen ≈ expected
end


@testset "centroid weights" begin
    x = [1. 4 10;
         0 3 0;
         -1 1 -11.3;]
    w = [1., 2, 0]

    expected = [3., 2., 1/3]

    cen = lguys.centroid(x, w)
    @test cen ≈ expected
end


@testset "centroid err" begin


end


@testset "centroid weights err" begin


end


@testset "centroid snap" begin
    N = 10
    pos = transpose(hcat(ones(N), zeros(N), -ones(N)))
    vel = transpose(hcat(zeros(N), 2*ones(N), zeros(N)))
    m = lguys.ConstVector(1., N)
    snap = lguys.Snapshot(pos=pos, vel=vel, m=m)

    cen = lguys.centroid(snap)

    @test cen.pos ≈ [1., 0., -1.]
    @test cen.vel ≈ [0., 2., 0.]
end




"""
given the number of points, the cdf of radii (f_x) and the cdf of velocities given r, f_v, creates a snapshot
"""
function make_snap(N, f_x, f_v, f_m=x->1.)
    rs = f_x.(rand(N))
    vs = f_v.(rs, rand(N))

    pos = lguys.rand_unit(N) * reshape(rs, 1, N)
    vel = lguys.rand_unit(N) * reshape(vs, 1, N)
    m = f_m.(rand(N))
    return lguys.Snapshot(pos=pos, vel=vel, m=m)
end


function uniform_snap(N=1000)
    f_x(p) = p
    f_v(r, p) = 0.1*p
    return make_snap(N, f_x, f_v)
end



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


@testset "uniform disk" begin

    f_x(p) = √p
    f_v(r, p) = 0.1*√p

end

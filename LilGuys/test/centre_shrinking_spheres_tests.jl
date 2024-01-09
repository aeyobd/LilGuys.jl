
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


function nfw_snap(N=1000)
    f_x(p) = p^2 / (1 + p)^2
    f_v(r, p) = 0.1 * p / (1 + p)
    return make_snap(N, f_x, f_v)
end


@testset "centre (ss)" begin
    for _ in 1:100
        cen = 20*randn(3)
        snap = uniform_snap()
        x_c, v_c, state = lguys.ss_centre(snap, cen)
        @test x_c ≈ cen
        @test v_c ≈ zeros(3)
    end
end

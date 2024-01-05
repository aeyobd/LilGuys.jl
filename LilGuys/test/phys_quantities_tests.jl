

function create_snapshot()
    pos = F[1  0  0;
            0 -1  0;
            0  0  1]

    vel = F[0 -2  0;
            1  0  0;
            0  0 -1.5]

    m = F[1, 1, 1]

    snap = lguys.Snapshot(pos=pos, vel=vel, m=m)
    return snap
end


@testset "r (matrix)" begin
    a = [1., 2., 3.]
    @test lguys.calc_r(a) ≈ 3.741657386773941

    b = [1. 0 1;
         2. 0 1;
         1. 0 -1]

    @test lguys.calc_r(b) ≈ [√6, 0, √3]
end


@testset "r (snap)" begin
    pos = F[1  0  0.5;
            0 -4  1.2;
            0  3  0.0]

    vel = randn(3, 3)

    m = F[1, 1, 1]

    snap = lguys.Snapshot(pos=pos, vel=vel, m=m)
    actual = lguys.calc_r.(snap)
    expected = [1, 5, 1.3]

    @test actual ≈ expected
end


@testset "kinetic energy" begin
    snap = create_snapshot()
    actual = lguys.calc_E_spec_kin(snap)
    expected = [0.5, 2, 1.125]

    @test actual ≈ expected
end


@testset "angular momentum" begin
    snap = create_snapshot()
    actual = lguys.calc_angular_momentum(snap)
    expected = F[0  0  0;
                 0  0  0;
                 1 -2  0]
                
    @test actual ≈ expected
    
end


@testset "E_spec" begin
    snap = create_snapshot()
end

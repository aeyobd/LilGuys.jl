using LilGuys


function create_snapshot()
    pos = F[1  0  0;
            0 -1  0;
            0  0  1]

    vel = F[0 -2  0;
            1  0  0;
            0  0 -1.5]

    m = F[1, 1, 1]

    snap = LilGuys.Snapshot(pos=pos, vel=vel, m=m)
    return snap
end


@testset "r (matrix)" begin
    a = [1., 2., 3.]
    @test LilGuys.r(a) ≈ 3.741657386773941

    b = [1. 0 1;
         2. 0 1;
         1. 0 -1]

    @test LilGuys.r(b) ≈ [√6, 0, √3]
end


@testset "r (snap)" begin
    pos = F[1  0  0.5;
            0 -4  1.2;
            0  3  0.0]

    vel = randn(3, 3)

    m = F[1, 1, 1]

    snap = LilGuys.Snapshot(pos=pos, vel=vel, m=m)
    actual = LilGuys.r.(snap)
    expected = F[1, 5, 1.3]

    @test actual ≈ expected
end


@testset "kinetic energy" begin
    snap = create_snapshot()
    actual = LilGuys.E_spec_kin(snap)
    expected = [0.5, 2, 1.125]

    @test actual ≈ expected
end


@testset "angular momentum" begin
    snap = create_snapshot()
    actual = LilGuys.angular_momentum(snap)
    expected = F[0  0  0;
                 0  0  0;
                 1 -2  0]
                
    println(actual)
    println(expected)
    @test actual ≈ expected
    
end


@testset "E_spec" begin
    snap = create_snapshot()
end


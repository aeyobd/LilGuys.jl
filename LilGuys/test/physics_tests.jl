function create_snapshot()
    pos = [1  0  0;
           0 -1  0;
           0  0  1]

    vel = [0 -2  0;
           1  0  0;
           0  0 -1.5]

    m = [1., 1, 1]

    snap = lguys.Snapshot(pos, vel, m)
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

@testset "r errors" begin
    a = [1., 2., 3., 4.]
    @test_throws DimensionMismatch lguys.calc_r(a)

    b = [1. 0 1 2;
         2. 0 1 2;
         1. 0 -1 4;
         1. 0 -1 5]

    @test_throws DimensionMismatch lguys.calc_r(b)
end


@testset "r (snap)" begin
    pos = [1.  0  0.5;
           0 -4  1.2;
           0  3  0.0]

    vel = randn(3, 3)

    m = [1., 1, 1]

    snap = lguys.Snapshot(pos, vel, m)
    actual = lguys.calc_r(snap)
    expected = [1, 5, 1.3]

    @test actual ≈ expected
end



@testset "kinetic energy" begin
    snap = create_snapshot()
    actual = lguys.calc_K_spec(snap)
    expected = [0.5, 2, 1.125]

    @test actual ≈ expected
end



@testset "angular momentum" begin
    snap = create_snapshot()
    actual = lguys.calc_L_spec(snap) .* snap.masses
    expected = [0.  0  0;
                0  0  0;
                1 -2  0]
                
    @test actual ≈ expected
    
    tot = lguys.calc_L_tot(snap)
    @test tot ≈ [0, 0, -1]
end



@testset "E_spec" begin
    Φs = [1., -0.23, 0.5, Inf]
    vs = [0.12, -0.2, π, 0]
end


@testset "circular v scalar" begin
    r = 0.
    M = Inf
    @test lguys.calc_v_circ(r, M) == 0.

    r = 1.
    M = 1
    @test lguys.calc_v_circ(r, M) == 1.

    r = 1.
    M = 2
    @test lguys.calc_v_circ(r, M) ≈ √2

    r = 2.
    M = 1 // 1
    @test lguys.calc_v_circ(r, M) ≈ 1/√2

    r = Inf
    M = Inf
    @test isnan(lguys.calc_v_circ(r, M))

    @test_throws DomainError lguys.calc_v_circ(1, -1)
    @test_throws DomainError lguys.calc_v_circ(-1, 1)
end


@testset "snap getters" begin
    N = 100
    pos = randn(3, N)
    vel = randn(3, N)
    m = rand(N)

    snap = lguys.Snapshot(pos, vel, m)

    @test lguys.get_x(snap) == pos[1, :]
    @test lguys.get_y(snap) == pos[2, :]
    @test lguys.get_z(snap) == pos[3, :]
    @test lguys.get_v_x(snap) == vel[1, :]
    @test lguys.get_v_y(snap) == vel[2, :]
    @test lguys.get_v_z(snap) == vel[3, :]
end



@testset "calc_W_tot" begin
    #  W only depends on masses and potentials

    N = 100
    pos = randn(3, N)
    m = rand(N)
    snap = Snapshot(pos, zeros(3, N), m)

    snap.Φs = lguys.calc_Φ(snap)

    W = 0

    for i in 1:N
        for j in i+1:N
            r = lguys.calc_r(pos[:, i], pos[:, j])
            W += m[i] * m[j] / r
        end
    end

    actual = lguys.calc_W_tot(snap)
    @test actual ≈ W
    @test actual ≈ -0.5 * sum(snap.Φs .* snap.masses)
end

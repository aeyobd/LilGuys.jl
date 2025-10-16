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
    @test lguys.radii(a) ≈ 3.741657386773941 

    b = [1. 0 1;
         2. 0 1;
         1. 0 -1]

    @test lguys.radii(b) ≈ [√6, 0, √3]
end


@testset "r (snap)" begin
    pos = [1.  0  0.5;
           0 -4  1.2;
           0  3  0.0]

    vel = randn(3, 3)

    m = [1., 1, 1]

    snap = lguys.Snapshot(pos, vel, m)
    actual = lguys.radii(snap)
    expected = [1, 5, 1.3]

    @test actual ≈ expected
end



@testset "kinetic energy" begin
    snap = create_snapshot()
    actual = lguys.kinetic_spec(snap)
    expected = [0.5, 2, 1.125]

    @test actual ≈ expected
end



@testset "angular momentum" begin
    snap = create_snapshot()
    actual = lguys.angular_momenta(snap) .* snap.masses
    expected = [0.  0  0;
                0  0  0;
                1 -2  0]
                
    @test actual ≈ expected
    
    tot = lguys.angular_momentum(snap)
    @test tot ≈ [0, 0, -1]
end



@testset "E_spec" begin
    Φs = [1., -0.23, 0.5, Inf]
    vs = [0.12, -0.2, π, 0]
end


@testset "circular v scalar" begin
    r = 0.
    M = Inf
    @test lguys.v_circ(r, M) == 0.

    r = 1.
    M = 1
    @test lguys.v_circ(r, M) == 1.

    r = 1.
    M = 2
    @test lguys.v_circ(r, M) ≈ √2

    r = 2.
    M = 1 // 1
    @test lguys.v_circ(r, M) ≈ 1/√2

    r = Inf
    M = Inf
    @test isnan(lguys.v_circ(r, M))

    @test_throws DomainError lguys.v_circ(1, -1)
    @test_throws DomainError lguys.v_circ(-1, 1)
end


@testset "snap getters" begin
    N = 100
    pos = randn(3, N)
    vel = randn(3, N)
    m = rand(N)

    snap = lguys.Snapshot(pos, vel, m)

    @test lguys.x_position(snap) == pos[1, :]
    @test lguys.y_position(snap) == pos[2, :]
    @test lguys.z_position(snap) == pos[3, :]
    @test lguys.x_velocity(snap) == vel[1, :]
    @test lguys.y_velocity(snap) == vel[2, :]
    @test lguys.z_velocity(snap) == vel[3, :]
end



@testset "W_tot" begin
    #  W only depends on masses and potentials

    N = 100
    pos = randn(3, N)
    m = rand(N)
    snap = Snapshot(pos, zeros(3, N), m)

    snap.potential = lguys.potential_nbody(snap)

    W = 0

    for i in 1:N
        for j in i+1:N
            r = lguys.radii(pos[:, i], pos[:, j])
            W += -m[i] * m[j] / r
        end
    end

    actual = lguys.potential_energy(snap)
    @test actual ≈ W
    @test actual ≈ 0.5 * sum(snap.potential .* snap.masses)
end

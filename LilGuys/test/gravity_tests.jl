@testset "f grav" begin
    positions = [0.;0;0;;]
    masses = ones(1)
    
    f(x) = lguys.calc_F_grav(positions, masses, x)

    x = [0.;0;0]
    @test f(x) ≈ [0.;0;0]

    x = [1.;0;0]
    @test f(x) ≈ [-1.;0;0]
end

@testset "Φ grav" begin
    pos = [0.;0;0;;]
    mass = ones(1)

    f(x_vec) = lguys.calc_Φ(pos, mass, x_vec)

    x_vec = [0.;0;0]
    @test f(x_vec) === -Inf

    x_vec = [1.;0;0]
    @test f(x_vec) ≈ -1.

    x_vec = [10.;0;0]
    @test f(x_vec) ≈ -1/10


    pos = [
           0. 2.
           0. 0.
           0. 0.
          ]
    mass = ones(2)
    f(x_vec) = lguys.calc_Φ(pos, mass, x_vec)
    @test f([1.;0;0]) ≈ -2.

    pos = [
           1. 0. 0. -4. -4. 0.
           0. 2. 0. 0. -3. 0.
           0. 0. 3. 0. 0. -6.
          ]

    mass = [1., 2, 3, 4, 5, 6]

    f(x_vec) = lguys.calc_Φ(pos, mass, x_vec)
    @test f([0.;0;0]) ≈ -6.

end

function make_rad_Φ(rs)
    println("here")
    N = length(rs)
    rs = reshape(rs, 1, N)
    positions = rs .* lguys.rand_unit(N) 
    masses = ones(N)
    velocities = zeros(3, N)
    snap = lguys.Snapshot(positions=positions, masses=masses, velocities=velocities)
    return lguys.calc_radial_Φ(snap)
end


@testset "radial Φ" begin
    f = make_rad_Φ([0.0, 1])
    @test f(0) ≈ -Inf
    @test f(1/2) ≈ -2 - 1
    @test f(1) ≈ -2
    @test f(2) ≈ -1
    @test f(10) ≈ -2/10

    f = make_rad_Φ([1, 2, 5])
    @test f(0) ≈ -1 - 1/2 - 1/5
    @test f(0.5) ≈ -1 - 1/2 - 1/5
    @test f(1) ≈ -1 - 1/2 - 1/5
    @test f(2) ≈ -1/2 - 1/2 - 1/5
    @test f(3) ≈ -1/3 - 1/3 - 1/5
    @test f(4) ≈ -1/4 - 1/4 - 1/5
    @test f(5) ≈ -3/5
    @test f(10) ≈ -3/10
end


@testset "radial Φ approx" begin
    N = 10000
    rs = rand(1, N) 
    positions = rs .* lguys.rand_unit(N) 
    masses = rand(N)
    velocities = zeros(3, N)
    snap =  lguys.Snapshot(positions=positions, masses=masses, velocities=velocities)

    rs = [0.0, 0.5, 1, 5, 10, 100]
    Φr = lguys.calc_radial_Φ(snap)
    f1(x) = Φr(lguys.calc_r(x))
    f(x) = lguys.calc_Φ(snap, x)

    rel_err(x) = abs(f(x) - f1(x)) / abs(f(x))

    @testset "rel err" begin
        for r in rs
            x = r * lguys.rand_unit()[:, 1]
            @test rel_err(x) < 3e-2
        end
    end
 
end


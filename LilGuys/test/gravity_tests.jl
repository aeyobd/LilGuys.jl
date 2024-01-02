@testset "f grav" begin
    pos = [0.;0;0;;]
    mass = ones(1)
    
    f(x) = lguys.calc_F_grav(pos, mass, x)

    x = [0.;0;0]
    @test f(x) ≈ [0.;0;0]

    x = [1.;0;0]
    @test f(x) ≈ [-1.;0;0]
end

@testset "Φ grav" begin
    pos = [0.;0;0;;]
    mass = ones(1)

    f(x) = lguys.calc_Φ(pos, mass, x)

    x = [0.;0;0]
    @test f(x) === -Inf

    x = [1.;0;0]
    @test f(x) ≈ -1.

    x = [10.;0;0]
    @test f(x) ≈ -1/10


    pos = [
           0. 2.
           0. 0.
           0. 0.
          ]
    mass = ones(2)
    f(x) = lguys.calc_Φ(pos, mass, x)
    @test f([1.;0;0]) ≈ -2.

    pos = [
           1. 0. 0. -4. -4. 0.
           0. 2. 0. 0. -3. 0.
           0. 0. 3. 0. 0. -6.
          ]

    mass = [1., 2, 3, 4, 5, 6]

    f(x) = lguys.calc_Φ(pos, mass, x)
    @test f([0.;0;0]) ≈ -6.

end

@testset "radial ϕ" begin
    N = 100000
    rs = rand(1, N) 
    pos = rs .* lguys.rand_unit(N) 
    m = rand(N)
    vel = zeros(3, N)
    snap =  lguys.Snapshot(pos=pos, vel=vel, m=m)

    f(x) = lguys.calc_Φ(pos, m, x) 

    rs = [0.0, 0.5, 1, 5, 10, 100]
    Φr = lguys.calc_radial_Φ(snap)
    f1(x) = Φr(lguys.calc_r(x))

    rel_err(x) = abs(f(x) - f1(x)) / abs(f(x))

    @testset "rel err" begin
        for r in rs
            x = r * lguys.rand_unit()[:, 1]
            println(r)
            @test rel_err(x) < 1e-3 
        end
    end
 
end


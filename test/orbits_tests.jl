tdir = mktempdir()


@testset "Orbit" begin
    @testset "construction"  begin
        pos = randn(3, 10)
        vel = randn(3, 10) / 2
        times = LinRange(2, 3.2, 10)

        orbit = Orbit(times=times, positions=pos, velocities=vel)
        @test LilGuys.positions(orbit) ≈ pos
        @test LilGuys.velocities(orbit) ≈ vel
        @test LilGuys.accelerations(orbit) |> isnothing

        accs = randn(3, 10)
        orbit = Orbit(times=times, positions=pos, velocities=vel, accelerations=accs)
        @test LilGuys.accelerations(orbit) ≈ accs

        @test orbit isa Orbit
    end

    @testset "improper construction" begin
        pos = randn(3, 9)
        vel = randn(3, 10)
        times = LinRange(3, 2.3, 10)
        @test_throws ArgumentError Orbit(positions=pos, velocities=vel, times=times) 

        pos = randn(3, 10)
        vel = randn(3, 10)
        times = LinRange(3, 2.3, 9)

        @test_throws ArgumentError Orbit(positions=pos, velocities=vel, times=times)

        pos = randn(3, 10)
        vel = randn(2, 10)
        times = LinRange(3, 2.3, 10)

        @test_throws ArgumentError Orbit(positions=pos, velocities=vel, times=times)

        pos = randn(4, 10)
        vel = randn(3, 10)
        times = LinRange(3, 2.3, 10)

        @test_throws  ArgumentError Orbit(positions=pos, velocities=vel, times=times)
    end


    @testset "orbit IO" begin
        orbit = Orbit(times=1:10, positions=randn(3, 10), velocities=randn(3, 10))
        write("$tdir/test.csv", orbit)
        orbit2 = Orbit("$tdir/test.csv")

        for attr in [LilGuys.positions, LilGuys.velocities, LilGuys.times, LilGuys.pericenter, LilGuys.apocenter, length]
            @test attr(orbit2) ≈ attr(orbit) atol=1e-8
        end

        @test isnothing(LilGuys.accelerations(orbit2))
        @test isnothing(LilGuys.accelerations(orbit))

        orbit = Orbit(times=1:10, positions=randn(3, 10), velocities=randn(3, 10), 
                      accelerations=rand(3, 10))

        write("$tdir/test.csv", orbit)
        orbit2 = Orbit("$tdir/test.csv")


        for attr in [LilGuys.positions, LilGuys.velocities, LilGuys.times, 
                     LilGuys.pericenter, LilGuys.apocenter, length, 
                     LilGuys.accelerations]
            @test attr(orbit2) ≈ attr(orbit) atol=1e-8
        end
    end


    @testset "to_frame" begin
        pos = randn(3, 8)
        vel = randn(3, 8) / 2
        times = LinRange(2, 3.2, 8)

        orbit = Orbit(times=times, positions=pos, velocities=vel)
        df = LilGuys.to_frame(orbit)
        @test df.t ≈ orbit.times
        @test df.x ≈ orbit.positions[1, :]
        @test df.y ≈ orbit.positions[2, :]
        @test df.z ≈ orbit.positions[3, :]
        @test df.v_x ≈ orbit.velocities[1, :]
        @test df.v_y ≈ orbit.velocities[2, :]
        @test df.v_z ≈ orbit.velocities[3, :]
        @test size(df) == (8, 7)

        accs = randn(3, 8)
        orbit = Orbit(times=times, positions=pos, velocities=vel, accelerations=accs)
        df = LilGuys.to_frame(orbit)

        @test df.t ≈ orbit.times
        @test df.x ≈ orbit.positions[1, :]
        @test df.y ≈ orbit.positions[2, :]
        @test df.z ≈ orbit.positions[3, :]
        @test df.v_x ≈ orbit.velocities[1, :]
        @test df.v_y ≈ orbit.velocities[2, :]
        @test df.v_z ≈ orbit.velocities[3, :]
        @test df.a_x ≈ orbit.accelerations[1, :]
        @test df.a_y ≈ orbit.accelerations[2, :]
        @test df.a_z ≈ orbit.accelerations[3, :]
    end
end

@testset "reverse" begin
    @testset "simple cases" begin
        pos = randn(3, 8)
        vel = randn(3, 8) / 2
        times = LinRange(2, 3.2, 8)

        orbit = Orbit(times=times, positions=pos, velocities=vel)
        orbit2 = reverse(orbit)

        @test reverse(orbit2.times) ≈ orbit.times
        @test orbit2.positions[:, end:-1:1] ≈ orbit.positions
        @test orbit2.velocities[:, end:-1:1] ≈ orbit.velocities
        @test length(orbit2) == 8
        @test LilGuys.pericenter(orbit2) ≈ LilGuys.pericenter(orbit)
        @test LilGuys.apocenter(orbit2) ≈ LilGuys.apocenter(orbit)


        accs = randn(3, 8)
        orbit = Orbit(times=times, positions=pos, velocities=vel, accelerations=accs)

        orbit2 = reverse(orbit)

        @test reverse(orbit2.times) ≈ orbit.times
        @test orbit2.positions[:, end:-1:1] ≈ orbit.positions
        @test orbit2.velocities[:, end:-1:1] ≈ orbit.velocities
        @test orbit2.accelerations[:, end:-1:1] ≈ orbit.accelerations
        @test length(orbit2) == 8
        @test LilGuys.pericenter(orbit2) ≈ LilGuys.pericenter(orbit)
        @test LilGuys.apocenter(orbit2) ≈ LilGuys.apocenter(orbit)


        pos = randn(3, 8)
        vel = randn(3, 8) / 2
        times = LinRange(2, 3.2, 8)

        orbit = Orbit(times=times, positions=pos, velocities=vel)
        orbit2 = reverse(orbit)
    end
end


@testset "leapfrog" begin
    @testset "isolation" begin
        f_acc(pos, vel, acc) = pos, vel, acc

    end

    @testset "kepler" begin

    end


    @testset "plummer" begin

    end

    @testset "disk" begin

    end
end


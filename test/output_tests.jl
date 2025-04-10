# this test requires quite a bit of setup
# as an output is a file of files :)
#

using HDF5 

let
    dir = mktempdir()
    mkdir(joinpath(dir, "out"))

    header = LilGuys.make_default_header(2, 0)
    snap1 = Snapshot(positions=[[0, 1, 2] [-1, 1, 0]], 
        velocities=[[0,0.1,0.2] [0.1,0.2,0.3]], masses=[0.5, 1], 
        index=[1,2], header=header)

    snap2 = Snapshot(positions=[[1,2,4] [2, 1, 2]], 
        velocities=[[0.1,0.1,0] [-0.1,0.1,0.1]], masses=[1, 0.5], 
        index=[2, 1], header=header)

    snap3 = Snapshot(positions=[[-2, 3.5, 0] [0., -1/2, 2]], 
        velocities=[[0., 0.15, 0.2] [0.05, 0.1, 0.3]], masses=[0.5, 1], 
        index=[1, 2], header=header)


    # write and open output
    LilGuys.write("$dir/out/snapshot_000.hdf5", snap1)
    LilGuys.write("$dir/out/snapshot_001.hdf5", snap2)
    LilGuys.write("$dir/out/snapshot_002.hdf5", snap3)

    old_dir = pwd()
    cd(dir * "/out")

    h5open("combined.hdf5", "w") do file
        HDF5.create_external(file, "/snap0", "snapshot_000.hdf5", "/")
        HDF5.create_external(file, "/snap1", "snapshot_001.hdf5", "/")
        HDF5.create_external(file, "/snap2", "snapshot_002.hdf5", "/")
    end

    cd(old_dir)
    out = Output(dir)


    @testset "output creation" begin
        @test length(out) == 3

        snaps = [snap1, snap2, snap3]
        for i in 1:3
            @test out[i].positions == snaps[i].positions
            @test out[i].velocities == snaps[i].velocities
            @test out[i].masses == snaps[i].masses
        end
    end


    @testset "snap int indexing" begin
        snap = Snapshot(joinpath(dir, "out/combined.hdf5/1"))
        @test snap.positions == snap1.positions
        @test snap.velocities == snap1.velocities
        @test snap.masses == snap1.masses
    end


    @testset "peris and apos" begin
        df = LilGuys.peris_apos(out)
        @test df.pericentre ≈ [√5, √2]
        @test df.apocentre ≈ [√16.25, √21]
    end


    @testset "extract" begin
        index = LilGuys.extract(out, :index)
        @test index == [[1, 2] [1, 2] [1, 2]]
    end

    @testset "extract with index" begin
        index = LilGuys.extract(out, :index, 1)
        @test index == [1,1,1]
    end

    @testset "extract vector" begin
        pos = LilGuys.extract_vector(out, :positions)

        @test pos[:, 1, 1] ≈ snap1.positions[:, 1]
        @test pos[:, 2, 1] ≈ snap1.positions[:, 2]
        @test pos[:, 1, 2] ≈ snap2.positions[:, 2]
        @test pos[:, 2, 2] ≈ snap2.positions[:, 1]
        @test pos[:, 1, 3] ≈ snap3.positions[:, 1]
        @test pos[:, 2, 3] ≈ snap3.positions[:, 2]

        vel1 = LilGuys.extract_vector(out, :velocities, 1)
        @test vel1[:, 1] ≈ snap1.velocities[:, 1]
        @test vel1[:, 2] ≈ snap2.velocities[:, 2]
    end


    @testset "print" begin
        s = sprint(print, out)
        @test s == "<output with 3 snapshots of 2 particles>"
    end
end


@testset "extract snapshot" begin
    pos  = [1 4
            2 5
            3 6]
    vel  = [0.1 0.2
            0.2 0.3
            0.3 0.4]
    mass = [0.5, 1]

    snap = Snapshot(pos, vel, mass, index=[2,1])

    @test LilGuys.extract_vector(snap, :positions) == pos[:, [2,1]]
    @test LilGuys.extract_vector(snap, :positions, 1) == pos[:, 2]

    @test LilGuys.extract(snap, :masses) == mass[[2,1]]
    @test LilGuys.extract(snap, :masses, 1) == mass[2]
end

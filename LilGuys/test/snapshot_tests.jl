using HDF5
using LilGuys


# Test for creating a default header
@testset "Snapshot Creation" begin
    N = 1000
    pos = randn((3, N))
    vel = randn((3, N))
    acc = randn((3, N))
    index = collect(1:N)
    Φ = -rand(N)
    Φ_ext = -rand(N)
    mass = 0
    h=rand()
    filename = "test.hdf5"
    header = LilGuys.make_default_header(N, mass)

    snap = Snapshot(pos=pos, vel=vel, acc=acc,
                    Φ=Φ, Φ_ext=Φ_ext, m=mass, index=index, h=h, 
                    filename=filename, header=header)


    @test all(snap.pos .== pos)
    @test all(snap.vel .== vel)
    @test all(snap.acc .== acc)
    @test all(snap.Φ .== Φ)
    @test snap.m == mass

    write!(filename, snap)

    snap1 = Snapshot(filename)
    @test all(snap1.pos .== snap.pos)
    @test all(snap1 .== snap)

end


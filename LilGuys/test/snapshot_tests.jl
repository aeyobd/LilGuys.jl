using HDF5
using LilGuys


function assert_equal(a::Snapshot, b::Snapshot)
    @test a.header == b.header
    @test a.m == b.m
    @test a.pos == b.pos
    @test a.vel == b.vel
    @test a.acc == b.acc
    @test a.Φ == b.Φ
end

# Test for creating a default header
@testset "Snapshot Creation" begin
    N = 1000
    pos = randn((3, N))
    vel = randn((3, N))
    acc = randn((3, N))
    index = collect(1:N)
    Φ = -rand(N)
    Φ_ext = -rand(N)
    m = LilGuys.ConstVector(0.1, N)
    h=rand()
    filename = "test.hdf5"
    header = LilGuys.make_default_header(N, m[1])

    snap = Snapshot(m=m, pos=pos, vel=vel,
                    acc=acc, Φ=Φ, Φ_ext=Φ_ext, index=index, h=h, 
                    filename=filename, header=header)


    @test all(snap.pos .== pos)
    @test all(snap.vel .== vel)
    @test all(snap.acc .== acc)
    @test all(snap.Φ .== Φ)
    @test snap.m[1] ≈ m[1]


    save(filename, snap)
    snap_saved = Snapshot(filename)
    assert_equal(snap, snap_saved)
end

function load_snap()
    return Snapshot("test.hdf5")
end

@testset "copy snapshot" begin
    snap = load_snap()
    snap2 = copy(snap)
    assert_equal(snap, snap2)

    N = 1000
    snap2.pos .+= [1,2,3]
    snap2.vel .-= [1,2,3]
    snap2.acc .*= 2
    snap2.Φ .*= 2
    println(typeof(snap2))
    println(typeof(snap))

    @test snap2.pos != snap.pos
    @test snap2.vel != snap.vel
    @test snap2.acc != snap.acc
    @test snap2.Φ != snap.Φ
end

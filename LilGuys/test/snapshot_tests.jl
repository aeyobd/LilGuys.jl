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
    snap1 = Snapshot(filename)
    for attr in fieldnames(Snapshot)
        if attr ∉ [:h, :header]
            @test all(getfield(snap, attr) .== getfield(snap1, attr))
        end
    end
    @test snap1.header == snap.header
end

function load_snap()
    return Snapshot("test.hdf5")
end

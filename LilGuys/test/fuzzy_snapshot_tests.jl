using HDF5
using LilGuys


function load_snap()
    return Snapshot("test.hdf5")
end

# Test for creating a default header
@testset "FSnapshot Creation" begin
    snap = load_snap()
    fsnap = LilGuys.FSnapshot(snap)

    @test fsnap.pos == snap.pos
    @test fsnap.vel == snap.vel
end




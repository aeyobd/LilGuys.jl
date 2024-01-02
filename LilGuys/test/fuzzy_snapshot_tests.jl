
function load_snap()
    return lguys.Snapshot("test.hdf5")
end

# Test for creating a default header
@testset "FSnapshot Creation" begin
    snap = load_snap()
    fsnap = lguys.FSnapshot(snap)

    @test fsnap.pos == snap.pos
    @test fsnap.vel == snap.vel
end






function assert_equal(a::Snapshot, b::Snapshot)
    for attr in fieldnames(Snapshot)
        if attr != :h
            @test getfield(a, attr) == getfield(b, attr)
        end
    end
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
    m = lguys.ConstVector(0.1, N)
    h=rand()
    filename = "test.hdf5"
    header = lguys.make_default_header(N, m[1])

    snap = Snapshot(masses=m, positions=pos, velocities=vel,
                    accelerations=acc, Φs=Φ, Φs_ext=Φ_ext, index=index, h=h, 
                    filename=filename, header=header)


    @test all(snap.positions .== pos)
    @test all(snap.velocities .== vel)
    @test all(snap.accelerations .== acc)
    @test all(snap.Φs .== Φ)
    @test snap.masses[1] ≈ m[1]


    lguys.save(filename, snap)
    snap_saved = Snapshot(filename)
    assert_equal(snap, snap_saved)
end


function load_snap()
    return Snapshot("test.hdf5")
end


@testset "copy snapshot" begin
    snap = load_snap()
    snap2 = deepcopy(snap)
    assert_equal(snap, snap2)

    N = 1000
    snap2.positions .+= [1,2,3]
    snap2.velocities .-= [1,2,3]
    snap2.accelerations .*= 2
    snap2.Φs .*= 2

    @test snap2.positions != snap.positions
    @test snap2.velocities != snap.velocities
    @test snap2.accelerations != snap.accelerations
    @test snap2.Φs != snap.Φs
end


# Test for creating a default header
@testset "Default Header Creation" begin
    N = 100
    mass = 1.0
    header = lguys.make_default_header(N, mass)

    @test header["NumPart_ThisFile"] == [0, N]
    @test header["NumPart_Total"] == [0, N]
    @test header["MassTable"] == [0.0, mass]

    expected_attrs = ["Time", "Redshift", "BoxSize", "NumFilesPerSnapshot"]
    @test all(k->k ∈ keys(header), expected_attrs)
end

# Test for getting and setting header attributes in an HDF5 file
@testset "HDF5 Header Get/Set" begin
    testfile = "test_header.h5"
    N = 100
    m = 1.0

    h5open(testfile, "w") do h5_f
        header = lguys.make_default_header(N, m)
        lguys.set_header!(h5_f, header)
    end

    h5open(testfile, "r") do h5_f
        header = lguys.get_header(h5_f)
        @test header["NumPart_ThisFile"] == [0, N]
        @test header["MassTable"] == [0.0, m]
    end

    rm(testfile)
end

# Test for getting and setting vectors in an HDF5 file
@testset "HDF5 Vector Get/Set" begin
    testfile = "test_vector.h5"
    vector_key = "Velocity"
    vector_val = [1.0, 2.0, 3.0]

    h5open(testfile, "w") do h5_f
        lguys.set_vector!(h5_f, vector_key, vector_val)
    end

    h5open(testfile, "r") do h5_f
        retrieved_vector = lguys.get_vector(h5_f, vector_key)
        @test all(retrieved_vector .== vector_val)
    end

    rm(testfile)
end


@testset "mass_is_fixed" begin
    N = 100
    pos = randn((3, N))
    vel = randn((3, N))
    masses = zeros(100)
    snap = Snapshot(pos, vel, masses)
    @test lguys.mass_is_fixed(snap) 

end


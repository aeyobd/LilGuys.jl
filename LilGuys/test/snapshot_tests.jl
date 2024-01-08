

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
    m = lguys.ConstVector(0.1, N)
    h=rand()
    filename = "test.hdf5"
    header = lguys.make_default_header(N, m[1])

    snap = Snapshot(m=m, pos=pos, vel=vel,
                    acc=acc, Φ=Φ, Φ_ext=Φ_ext, index=index, h=h, 
                    filename=filename, header=header)


    @test all(snap.pos .== pos)
    @test all(snap.vel .== vel)
    @test all(snap.acc .== acc)
    @test all(snap.Φ .== Φ)
    @test snap.m[1] ≈ m[1]


    lguys.save(filename, snap)
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

    @test snap2.pos != snap.pos
    @test snap2.vel != snap.vel
    @test snap2.acc != snap.acc
    @test snap2.Φ != snap.Φ
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


# Test for the epsilon value extraction
@testset "Epsilon Value Extraction" begin
    # Create a test file with the expected content
    ϵ = 0.0512
    directory = mktempdir()
    testfile = joinpath(directory, "parameters-usedvalues")
    open(testfile, "w") do file
        write(file, "SomeParameter 1.0\n")
        write(file, "SofteningComovingClass0 $(ϵ)\n")
        write(file, "AnotherParameter 2.0\n")
    end

    # Test the epsilon value extraction function
    epsilon = lguys.get_epsilon(directory)
    @test epsilon === ϵ

    # Cleanup the temporary directory
    rm(directory, recursive=true)
end


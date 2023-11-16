using HDF5
using LilGuys

# Test for creating a default header
@testset "Default Header Creation" begin
    N = 100
    mass = 1.0
    header = LilGuys.make_default_header(N, mass)

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
        header = LilGuys.make_default_header(N, m)
        LilGuys.set_header!(h5_f, header)
    end

    h5open(testfile, "r") do h5_f
        header = LilGuys.get_header(h5_f)
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
        LilGuys.set_vector!(h5_f, vector_key, vector_val)
    end

    h5open(testfile, "r") do h5_f
        retrieved_vector = LilGuys.get_vector(h5_f, vector_key)
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
    epsilon = LilGuys.get_epsilon(directory)
    @test epsilon === ϵ

    # Cleanup the temporary directory
    rm(directory, recursive=true)
end



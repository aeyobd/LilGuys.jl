
tdir = mktempdir()

function assert_equal(a::Snapshot, b::Snapshot)
    for attr in fieldnames(Snapshot)
        if attr != :h
            if getfield(a, attr) === NaN
                @test isnan(getfield(b, attr))
            else
                @test getfield(a, attr) == getfield(b, attr)
            end
        end
    end
end


function assert_approx(a::Snapshot, b::Snapshot)
    for attr in fieldnames(Snapshot)
        if attr != :h
            if getfield(a, attr) === NaN
                @test isnan(getfield(b, attr))
            else
                x = getfield(a, attr)
                y = getfield(b, attr)

                if x isa AbstractArray
                    @test x ≈ y
                elseif x isa Dict
                    # pass
                else
                    @test getfield(a, attr) == getfield(b, attr)
                end
            end
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
    filename = joinpath(tdir, "test.hdf5")
    header = lguys.make_default_header(N, m[1])

    snap = Snapshot(masses=m, positions=pos, velocities=vel,
                    accelerations=acc, potential=Φ, potential_ext=Φ_ext, index=index, h=h, 
                    filename=filename, header=header)


    @test all(snap.positions .== pos)
    @test all(snap.velocities .== vel)
    @test all(snap.accelerations .== acc)
    @test all(snap.potential .== Φ)
    @test all(snap.potential_ext .== Φ_ext)
    @test snap.masses[1] ≈ m[1]
    @test size(snap) == (N,)


    lguys.write(filename, snap)
    snap_saved = Snapshot(filename)
    assert_equal(snap, snap_saved)
end


@testset "simple snapshot creation" begin
    N = 100
    pos = randn((3, N))
    vel = randn((3, N))
    mass = 2.5
    snap = Snapshot(pos, vel, mass)

    @test all(snap.positions .== pos)
    @test all(snap.velocities .== vel)
    @test all(snap.masses .== mass)

    masses = rand(N)
    snap = Snapshot(pos, vel, masses)
    @test all(snap.masses .== masses)

    @test_throws DimensionMismatch Snapshot(pos, vel, rand(N+1))
    @test_throws DimensionMismatch Snapshot(pos, vel[:, 1:N-2], masses)
end



function load_snap()
    return Snapshot(joinpath(tdir, "test.hdf5"))
end


@testset "copy snapshot" begin
    snap = load_snap()
    snap2 = deepcopy(snap)
    assert_equal(snap, snap2)

    N = 1000
    snap2.positions .+= [1,2,3]
    snap2.velocities .-= [1,2,3]
    snap2.accelerations .*= 2
    snap2.potential .*= 2

    @test snap2.positions != snap.positions
    @test snap2.velocities != snap.velocities
    @test snap2.accelerations != snap.accelerations
    @test snap2.potential != snap.potential
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
    testfile = joinpath(tdir, "test_header.hdf5")
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
    testfile = joinpath(tdir, "test_vector.hdf5")
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

@testset "regenerate_header!" begin
    # mass and particle number stored in header, want to 
    # make sure this is updated if we save a new snapshot
    N = 4
    snap = Snapshot(randn((3, N)), randn((3, N)), 1.5)
    snap.masses = [1,2,3,4]
    
    lguys.regenerate_header!(snap)
    lguys.write(joinpath(tdir, "test.hdf5"), snap)
    snap2 = Snapshot(joinpath(tdir, "test.hdf5"))
    @test snap2.masses == snap.masses


    snap.masses = ones(4)
    
    #lguys.regenerate_header!(snap)
    # check this is done automatically
    lguys.write(joinpath(tdir, "test.hdf5"), snap)
    snap2 = Snapshot(joinpath(tdir, "test.hdf5"))
    @test snap2.masses == lguys.ConstVector(1.0, 4)

    snap = snap[1:3]
    #lguys.regenerate_header!(snap)
    lguys.write(joinpath(tdir, "test.hdf5"), snap)
    snap2 = Snapshot(joinpath(tdir, "test.hdf5"))
    @test snap2.masses == lguys.ConstVector(1.0, 3)
    @test snap2.header["NumPart_ThisFile"] == [0, 3]
end

@testset "rescale" begin
    N = 100
    pos = randn((3, N))
    vel = randn((3, N))
    masses = ones(N)
    snap = Snapshot(pos, vel, masses)

    m_scale = 2
    r_scale = 0.1
    snap = lguys.rescale(snap, m_scale, r_scale)
    v_scale = sqrt(m_scale / r_scale)

    snap_rescaled = Snapshot(pos*r_scale, vel*v_scale, masses*m_scale)

    assert_approx(snap, snap_rescaled)

    # TODO: add potential tests....

end



@testset "add_stars!" begin
    pos = [1 2 3
           2 3 4
           1 0 1]
    vel = zeros(3, 3)

    idx = [3,1,2]
    snap = Snapshot(pos, vel, 1.0, index=[4,1,2])

    probabilities = [0.2, 0.1, 0.3]
    lguys.add_stars!(snap, probabilities)
    @test snap.weights ≈ [0.3, 0.2, 0.1]

    lguys.add_stars!(snap, [1,2,4], probabilities)
    @test snap.weights ≈ [0.3, 0.2, 0.1]

    @test_throws DimensionMismatch lguys.add_stars!(snap, [1, 0.2])
    @test_throws AssertionError lguys.add_stars!(snap, [1,2,3], [0.1, 0.2])
    @test_throws AssertionError lguys.add_stars!(snap, [1,2], [0.1, 0.2])
end

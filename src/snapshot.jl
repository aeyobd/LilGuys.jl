using Printf
import Base: @kwdef


"""A map of the HDF5 columns to the snapshot fields"""
const h5vectors = Dict(
    :Φs=>"Potential",
    :Φs_ext=>"ExtPotential",
    :index=>"ParticleIDs",
    :positions=>"Coordinates",
    :velocities=>"Velocities",
    :accelerations=>"Acceleration",
    :masses=>"Masses",
   )

const snap_matricies = [:positions, :velocities, :accelerations]
const snap_vectors = [:masses, :Φs, :Φs_ext, :index]



"""
A snapshot of a gadget simulation. Units are all code units.
"""
@kwdef mutable struct Snapshot 
    """The positions of the particles"""
    positions::Matrix{F}

    """The velocities of the particles"""
    velocities::Matrix{F}

    """The masses of the particles"""
    masses::Union{Vector, ConstVector}

    """The particle IDs"""
    index::Vector{Int}  

    """The accelerations of the particles"""
    accelerations::OptMatrix = nothing

    """The potential of the particles"""
    Φs::OptVector = nothing  

    """The external (milky way) potential of the particles"""
    Φs_ext::OptVector = nothing

    """The filename of the snapshot"""
    filename::String = ""

    """softening length"""
    h::Real = NaN

    """The Gadget HDF5 header"""
    header::Dict{String, Any}

    """The (adopted) centre"""
    x_cen::Vector{F} = zeros(F, 3)

    """The (adopted) velocity centre"""
    v_cen::Vector{F} = zeros(F, 3)

    """The (stellar) weights of the particles"""
    weights::Union{Vector, ConstVector} = ConstVector(1.0, length(index))

    time::F = NaN

    """The radii of the particles, stored on first calculation"""
    _radii::OptVector = nothing

end




"""
    Snapshot(positions, velocities, mass::Real)

Create a snapshot with constant mass.
"""
function Snapshot(positions, velocities, mass::Real; kwargs...)
    N = size(positions, 2)
    masses = ConstVector(mass, N)
    return Snapshot(positions, velocities, masses; kwargs...)
end


"""
    Snapshot(positions, velocities, masses)

Create a snapshot.
"""
function Snapshot(positions, velocities, masses; time=0,kwargs...)
    N = size(positions, 2)
    if size(velocities) != size(positions)
        throw(DimensionMismatch("velocities and positions must have the same size"))
    end

    if mass_is_fixed(masses)
        m_header = masses[1]
    else
        m_header = 0.0
        if length(masses) != N
            throw(DimensionMismatch("masses must have the same length as the number of particles"))
        end
    end

    header = make_default_header(N, m_header, time)
    index = collect(1:N)
    return Snapshot(positions=positions, velocities=velocities, masses=masses, header=header, index=index; time=time,kwargs...)
end



"""
    Snapshot(filename)

Load a snapshot from an HDF5 file.
The filename may be a snapshot.hdf5 file or can be the path to an output
with a slash-index to retrieve the ith snapshot from the output (including the
centre).
"""
function Snapshot(filename::String)
    outidx_pattern = r"/(-?\d+)$"

    local snap

    if splitext(filename)[2] == ".hdf5"
        h5open(filename, "r") do h5f
            snap = Snapshot(h5f, filename=filename)
        end
    elseif occursin(outidx_pattern, filename)
        base_path = replace(filename, outidx_pattern => "")
        num = parse(Int, match(outidx_pattern, filename).captures[1])
        out = Output(base_path)
        if num < 0
            num = length(out) + 1 + num
        end
        @info "loading snapshot $(out.index[num]) from $base_path"
        snap = out[num]
    else
        throw(ArgumentError("filename must be an HDF5 file or an output/index path"))
    end

    return snap
end


"""
    Snapshot(h5f; mmap=false, filename)

Load a snapshot from an HDF5 file.
"""
function Snapshot(h5f::HDF5.H5DataStore; mmap=false, filename="", group="PartType1")
    kwargs = Dict{Symbol, Any}()

    for (var, col) in h5vectors
        if col ∈ keys(h5f[group])
            kwargs[var] = get_vector(h5f, "$group/$col", mmap=mmap)
        end
    end

    header = get_header(h5f)
    kwargs[:header] = header
    m = header["MassTable"][2]
    N = header["NumPart_ThisFile"][2]

    if m != 0  || (m == 0 && :masses ∉ keys(kwargs))
        kwargs[:masses] = ConstVector(m, N)
    end

    kwargs[:filename] = filename

    return Snapshot(; kwargs...)
end


function Base.size(snap::Snapshot) 
    return (length(snap.index),)
end


function Base.length(snap::Snapshot) 
    return length(snap.index)
end

# Base.IndexStyle(::Type{<:Snapshot}) = IndexLinear()


function Base.getindex(snap::Snapshot, idx)
    kwargs = Dict{Symbol, Any}()
    kwargs[:h] = snap.h
    kwargs[:x_cen] = snap.x_cen
    kwargs[:v_cen] = snap.v_cen
    kwargs[:header] = snap.header
    kwargs[:filename] = snap.filename
    kwargs[:time] = snap.time

    for sym in [:positions, :velocities, :masses, :index, :accelerations, :Φs, :Φs_ext, :weights]
        if getproperty(snap, sym) === nothing
            continue
        elseif sym ∈ snap_matricies
            kwargs[sym] = getproperty(snap, sym)[:, idx]
        elseif sym ∈ snap_vectors
            kwargs[sym] = getproperty(snap, sym)[idx]
        elseif sym ∈ [:weights]
            if isa(getproperty(snap, sym), ConstVector)
                kwargs[sym] = getproperty(snap, sym)
            else
                kwargs[sym] = getproperty(snap, sym)[idx]
            end
        else
            continue
        end
    end

    return Snapshot(; kwargs...)
end



"""is the particle mass constant in the snapshot?"""
function mass_is_fixed(snap::Snapshot)
    return mass_is_fixed(snap.masses)
end


function mass_is_fixed(masses::Union{Vector, ConstVector})
    return masses isa ConstVector && masses[1] != 0 || all(masses .== masses[1])
end


"""
    save(filename, snap)

Save a snapshot to an HDF5 file.

Notes: does regenerate the snapshot's header.
"""
function save(filename::String, snap::Snapshot; centre=false)
    h5open(filename, "w") do h5f
        if centre
            snap1 = deepcopy(snap)
            snap1.positions .-= snap1.x_cen
            snap1.velocities .-= snap1.v_cen
        else
            snap1 = snap
        end

        save(h5f, snap1)
    end
end


function save(snap::Snapshot)
    save(snap.filename, snap)
end



function save(h5f::HDF5.H5DataStore, snap::Snapshot)
    regenerate_header!(snap)

    set_header!(h5f, snap.header)

    for (var, _) in h5vectors
        save_vector(h5f, snap, var)
    end
end


"""
saves a vector from a snapshot into an HDF5 file
"""
function save_vector(h5f::HDF5.H5DataStore, snap::Snapshot, var::Symbol; group="PartType1")
    col = h5vectors[var]
    val = getproperty(snap, var) 

    if group ∉ keys(h5f)
        HDF5.create_group(h5f, group)
    end

    if var == :masses && mass_is_fixed(snap)
        # pass
    elseif val !== nothing
        set_vector!(h5f, "$group/$col", val)
    end
end



function Base.show(io::IO, snap::Snapshot)
    print(io, "<snapshot with $(length(snap)) particles>")
    return io
end



"""
    make_default_header(N, mass)

Create a default Gadget header for a snapshot with N particles and DM mass `mass`.
"""
function make_default_header(N, mass, time=0)

    return Dict(
        "NumPart_ThisFile"=>UInt64[0, N],
        "MassTable"=>F[0.0, mass],
        "Time"=>time,
        "Redshift"=>0.0,
        "NumPart_Total"=>UInt64[0, N],
        "NumFilesPerSnapshot"=>1,
        "BoxSize"=>1000,
       )

end


"""
Regenerates the header of the snapshot, accounting for changes in particle mass and number of particles.
"""
function regenerate_header!(snap::Snapshot)
    if mass_is_fixed(snap)
        snap.header["MassTable"][2] = snap.masses[1]
    else
        snap.header["MassTable"][2] = 0.0
    end

    if length(snap) != snap.header["NumPart_ThisFile"][2]
        snap.header["NumPart_ThisFile"][2] = length(snap)
    end
end


"""
    make_gadget2_header(N, mass)

Create a Gadget-2 header for a snapshot with N particles and DM mass `mass`.
"""
function make_gadget2_header(N, mass)
    return Dict(
        "NumPart_ThisFile"=>F[0, N, 0, 0, 0, 0],
        "MassTable"=>F[0.0, mass, 0.0, 0.0, 0.0, 0.0],
        "Time"=>0.0,
        "Redshift"=>0.0,
        "NumPart_Total"=>F[0, N, 0, 0, 0, 0],
        "NumFilesPerSnapshot"=>1,
        "Flag_Entropy_ICs"=>0,
        "NumPart_Total_HighWord"=>F[0, 0, 0, 0, 0, 0],
    )
end


"""
    rescale(snapshot, mass_scale, radius_scale)

Returns a snapshot scaled by the given mass and radius (using code units)
"""
function rescale(snap::Snapshot, m_scale::Real, r_scale::Real)
    v_scale = sqrt(G * m_scale / r_scale)
    a_scale = G*m_scale / r_scale^2 # = v_scale^2 / r_scale
    Φ_scale = G * m_scale / r_scale # = v_scale^2

    rescaled = deepcopy(snap)

    rescaled.positions = snap.positions * r_scale
    rescaled.velocities = snap.velocities * v_scale
    rescaled.masses = snap.masses * m_scale

    if snap.Φs !== nothing
        rescaled.Φs .*= Φ_scale
    end

    if snap.Φs_ext !== nothing
        rescaled.Φ_ext .*= Φ_scale
    end


    if snap.accelerations !== nothing
        rescaled.accelerations .*= a_scale
    end

    return rescaled
end


"""
    add_stars!(snap, [index, ]probability)

Given a snapshot, add a set of stars with the given probability.
Probabilites should be sorted in the index order of the snapshot.
If index is given, this is interpreted as the indices of the stars of the
probability array, and the snapshot index is checked against the provided index.
"""
function add_stars!(snap::Snapshot, index, probability)
    @assert sort(snap.index) == index "stars must have the same indices as the snapshot"

    add_stars!(snap, probability)
    return snap
end


function add_stars!(snap::Snapshot, probability)
    if length(probability) != length(snap)
        throw(DimensionMismatch("stars must have the same length as the snapshot"))
    end

    if isperm(snap.index)
        idx = snap.index
    else
        idx = invperm(sortperm(snap.index))
    end

    snap.weights = probability[idx]
    return snap
end

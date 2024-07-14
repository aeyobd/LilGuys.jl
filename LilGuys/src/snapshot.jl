using LinearAlgebra: norm
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
    weights::Union{Vector, ConstVector} = ConstVector(1.0, 0)

    """The radii of the particles, stored on first calculation"""
    _radii::OptVector = nothing
end




"""
    Snapshot(positions, velocities, mass::Real)

Create a snapshot with constant mass.
"""
function Snapshot(positions, velocities, mass::Real)
    N = size(positions, 2)
    masses = ConstVector(mass, N)
    return Snapshot(positions, velocities, masses)
end


"""
    Snapshot(positions, velocities, masses)

Create a snapshot.
"""
function Snapshot(positions, velocities, masses)
    N = size(positions, 2)
    if mass_is_fixed(masses)
        m_header = masses[1]
    else
        m_header = 0.0
    end

    header = make_default_header(N, m_header)
    index = collect(1:N)
    return Snapshot(positions=positions, velocities=velocities, masses=masses, header=header, index=index)
end



"""
    Snapshot(filename)

Load a snapshot from an HDF5 file.
"""
function Snapshot(filename::String)
    h5open(filename, "r") do h5f
        return Snapshot(h5f, filename=filename)
    end
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
    for sym in fieldnames(Snapshot)
        if getproperty(snap, sym) === nothing
            continue
        elseif !isa(idx, Int) && sym ∈ [:header, :filename]
            kwargs[sym] = getproperty(snap, sym)
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
"""
function save(filename::String, snap::Snapshot)
    h5open(filename, "w") do h5f
        save(h5f, snap)
    end
end


function save(snap::Snapshot)
    save(snap.filename, snap)
end



function save(h5f::HDF5.H5DataStore, snap::Snapshot)
    if mass_is_fixed(snap)
        snap.header["MassTable"][2] = snap.masses[1]
    end
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
function make_default_header(N, mass)

    return Dict(
        "NumPart_ThisFile"=>F[0, N],
        "MassTable"=>F[0.0, mass],
        "Time"=>0.0,
        "Redshift"=>0.0,
        "NumPart_Total"=>F[0, N],
        "NumFilesPerSnapshot"=>1,
        "BoxSize"=>1000,
       )

end


"""
Regenerates the header of the snapshot.
"""
function regenerate_header!(snap::Snapshot)
    if mass_is_fixed(snap)
        m = snap.masses[1]
    else
        m = 0.0
    end

    N = length(snap)
    snap.header = make_default_header(N, m)
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



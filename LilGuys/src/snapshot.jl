using HDF5
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

Attributes
----------
positions : 3xN matrix
    The positions of the particles
velocities : 3xN matrix
    The velocities of the particles
masses : N vector
    The masses of the particles
index : N vector of Int
    The particle IDs

accelerations : 3xN matrix (optional)
    The accelerations of the particles
Φs : N vector (optional)
    The potential of the particles
Φs_ext : N vector (optional)
    The external potential of the particles

filename : str
    The filename of the snapshot

h : Real
    softening length
header : Dict{String, Any}
    The Gadget HDF5 header
x_cen : 3 vector
    The (adopted) centre
v_cen : 3 vector
    The (adopted) velocity centre


radii : 3 vector (optional)
    The radii of the particles, stored on first calculation
"""
@kwdef mutable struct Snapshot 
    positions::Matrix{F}
    velocities::Matrix{F}
    masses::Union{Vector, ConstVector}
    index::Vector{Int}  

    accelerations::OptMatrix = nothing
    Φs::OptVector = nothing  
    Φs_ext::OptVector = nothing

    filename::String = ""
    h::Real = NaN
    header::Dict{String, Any}

    x_cen::Vector{F} = zeros(F, 3)
    v_cen::Vector{F} = zeros(F, 3)
    weights::Union{Vector, ConstVector} = ConstVector(1.0, 0)

    # cahced
    #
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



function Snapshot(filename::String)
    h5open(filename, "r") do h5f
        return Snapshot(h5f, filename=filename)
    end
end


function Snapshot(h5f::HDF5.H5DataStore; mmap=false, filename="")
    kwargs = Dict{Symbol, Any}()

    for (var, col) in h5vectors
        if col ∈ keys(h5f["PartType1"])
            kwargs[var] = get_vector(h5f, col, mmap=mmap)
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


Base.size(snap::Snapshot) = (length(snap.index),)
Base.length(snap::Snapshot) = length(snap.index)
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


function Base.copy(snap::Snapshot)
    kwargs = Dict{Symbol, Any}()
    for sym in fieldnames(Snapshot)
        kwargs[sym] = deepcopy(getproperty(snap, sym))
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




function Base.show(io::IO, snap::Snapshot)
    print(io, "<snapshot with $(length(snap)) particles>")
    return io
end


##### HDF5 functions #####


# Make a default header for an HDF5 file
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

function regenerate_header!(snap::Snapshot)
    if mass_is_fixed(snap)
        m = snap.masses[1]
    else
        m = 0.0
    end

    N = length(snap)
    snap.header = make_default_header(N, m)
end


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


# Set the header in an HDF5 file
function set_header!(h5f::HDF5.File, header::Dict{String,Any})
    if "Header" ∉ keys(h5f)
        create_group(h5f, "Header")
    end
    h5_header = h5f["Header"]
    for (key, val) in header
        set_header_attr(h5f, key, val)
    end
end


# gets the gadget header of an HDF5 file
function get_header(h5f::HDF5.H5DataStore)
    return Dict(attrs(h5f["Header"]))
end


# gets a vector from an HDF5 file
function get_vector(h5f::HDF5.H5DataStore, key::String; mmap=false, group="PartType1")
    path = group * "/" * key
    if mmap
        return HDF5.readmmap(h5f[path])
    else
        return read(h5f[path])
    end
end


function set_vector!(h5f::HDF5.File, key::String, val, group="PartType1")
    if group ∉ keys(h5f)
        create_group(h5f, group)
    end
    path = group * "/" * key
    h5f[path] = val
end


function set_vector_ele!(h5f::HDF5.File, key::String, el::Int, val, group="PartType1")
    path = group * "/" * key
    h5f[path][el] = val
end


function set_header_attr(h5f::HDF5.File, key::String, val)
    header = attrs(h5f["Header"])
    header[key] = val
end


function save_vector(h5f::HDF5.H5DataStore, snap::Snapshot, var::Symbol)
    col = h5vectors[var]
    val = getproperty(snap, var) 

    if var == :masses && mass_is_fixed(snap)
        # pass
    elseif val !== nothing
        set_vector!(h5f, col, val)
    end
end


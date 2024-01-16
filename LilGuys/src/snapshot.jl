using HDF5
using LinearAlgebra: norm
using Printf
import Base: @kwdef




# this maps the gadget vectors to the snapshot attributes
const h5vectors = Dict(
    :Φs=>"Potential",
    :Φs_ext=>"ExtPotential",
    :index=>"ParticleIDs",
    :positions=>"Coordinates",
    :velocities=>"Velocities",
    :accelerations=>"Acceleration",
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

"""
@kwdef struct Snapshot 
    positions::Matrix{F}
    velocities::Matrix{F}
    masses::Union{Vector{F}, ConstVector}
    accelerations::OptMatrix = nothing
    Φs::OptVector = nothing  # potentials
    Φs_ext::OptVector = nothing
    header::Dict{String, Any} = make_default_header(size(positions, 2), masses[1]) #TODO: more robustly deal with variable masses
    index::Vector{Int} = collect(1:size(positions, 2))
    filename::String = ""
    h::Real = NaN
end


function Snapshot(filename::String; mmap=false)
    kwargs = Dict{Symbol, Any}()

    h5open(filename, "r") do h5f
        for (var, col) in h5vectors
            if col ∈ keys(h5f["PartType1"])
                kwargs[var] = get_vector(h5f, col, mmap=mmap)
            end

        end
        header = get_header(h5f)
        kwargs[:header] = header
        m = header["MassTable"][2]
        kwargs[:masses] = ConstVector(m, size(kwargs[:positions], 2))
        kwargs[:filename] = filename
    end
    return Snapshot(; kwargs...)
end


Base.size(snap::Snapshot) = (length(snap.index),)
Base.length(snap::Snapshot) = length(snap.index)
# Base.IndexStyle(::Type{<:Snapshot}) = IndexLinear()



function Base.getindex(snap::Snapshot, idx::Union{UnitRange, Vector, Colon, BitVector, Int})
    kwargs = Dict{Symbol, Any}()
    kwargs[:h] = snap.h
    for sym in fieldnames(Snapshot)
        if getproperty(snap, sym) === nothing
            continue
        elseif !isa(idx, Int) && sym ∈ [:header, :filename]
            kwargs[sym] = getproperty(snap, sym)
        elseif sym ∈ snap_matricies
            kwargs[sym] = getproperty(snap, sym)[:, idx]
        elseif sym ∈ snap_vectors
            kwargs[sym] = getproperty(snap, sym)[idx]
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


function save(snap::Snapshot)
    save(snap.filename, snap)
end


function save(filename::String, snap::Snapshot)
    h5open(filename, "w") do h5f
        for (var, col) in h5vectors
            val = getproperty(snap, var) # matrix is unneccesary,
            if val === nothing
                println("warning: $var is nothing")
                if col ∈ ["Potential", "ExtPotential"]
                    val = zeros(F, length(snap))
                else
                    val = zeros(F, 3, length(snap))
                end
            end
            set_vector!(h5f, col, val)
        end

        set_header!(h5f, snap.header)
    end
end


function Base.show(io::IO, snap::Snapshot)
    print(io, "<snapshot with $(length(snap)) particles>")
    return io
end


##### HDF5 functions #####


# Make a default header for an HDF5 file
function make_default_header(N, mass)
    header = Dict{String,Any}()
    header["NumPart_ThisFile"] = F[0, N]
    header["NumPart_Total"] = F[0, N]
    header["MassTable"] = F[0.0, mass]
    header["Time"] = 0.0
    header["Redshift"] = 0.0
    header["BoxSize"] = 350.0
    header["NumFilesPerSnapshot"] = 1
    return header
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
function get_header(h5f::HDF5.File)
    return Dict(attrs(h5f["Header"]))
end

# gets a vector from an HDF5 file
function get_vector(h5f::HDF5.File, key::String; mmap=false, group="PartType1")
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



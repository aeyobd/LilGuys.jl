using HDF5
using LinearAlgebra: norm
using Printf


# this maps the gadget vectors to the snapshot attributes
const h5vectors = Dict(
    :Φ=>"Potential",
    :Φ_ext=>"ExtPotential",
    :index=>"ParticleIDs",
    :pos=>"Coordinates",
    :vel=>"Velocities",
    :acc=>"Acceleration",
   )


AbstractSnapshot = AbstractArray{Particle, 1}

Base.@kwdef struct Snapshot <: AbstractSnapshot
    pos::Matrix{F}
    vel::Matrix{F}
    m::Vector{F}
    acc::Matrix{F} = []
    Φ::Vector{F} = []
    Φ_ext::Vector{F} = []
    header::Dict{String, Any} = make_default_header(size(pos, 2), m[1]) #TODO: more robustly deal with variable masses
    index::Vector{Int} = collect(1:size(pos, 2))
    filename::String = ""
    h::F = NaN
end


function Snapshot(filename::String; mmap=false)
    kwargs = Dict{Symbol, Any}()

    h5open(filename, "r") do h5f
        for (var, header) in h5vectors
            if header ∈ keys(h5f["PartType1"])
                kwargs[var] = get_vector(h5f, header, mmap=mmap)
            end
        end
        header = get_header(h5f)
        kwargs[:header] = header
        m = header["MassTable"][2]
        kwargs[:m] = ConstVector(m, size(kwargs[:pos], 2))
        kwargs[:filename] = filename
    end
    return Snapshot(; kwargs...)
end


Base.size(snap::Snapshot) = (length(snap.index),)
Base.IndexStyle(::Type{<:Snapshot}) = IndexLinear()


function iloc(snap::Snapshot, i::Int)
    # return sortperm(snap.index)[i]
end


function Base.getindex(snap::Snapshot, idx::Union{UnitRange, Vector, Colon, BitVector, Int})
    kwargs = Dict{Symbol, Any}()
    kwargs[:h] = snap.h
    for sym in [:m, :Φ, :Φ_ext, :index]
        kwargs[sym] = getproperty(snap, sym)[idx]
    end
    for sym in [:pos, :vel, :acc]
        kwargs[sym] = getproperty(snap, sym)[:, idx]
    end

    if isa(idx, Int)
        return Particle(; kwargs...)
    end

    kwargs[:header] = snap.header
    kwargs[:filename] = snap.filename
    return Snapshot(; kwargs...)
end


function Base.setindex!(snap::Snapshot, p::Particle, i::Int)
    kwargs = Dict{Symbol, Any}()

    for sym in [:m, :pos, :vel, :Φ, :Φ_ext, :acc]
        getproperty(snap, sym)[i] = getproperty(particle, sym)
    end
    return snap[idx]
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
            set_vector!(h5f, col, val)
        end

        set_header!(h5f, snap.header)
    end
end


function Base.show(io::IO, snap::Snapshot)
    print(io, "<snapshot with $(length(snap)) particles>")
    return io
end



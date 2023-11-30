using HDF5
using LinearAlgebra: norm
using Printf


const h5vectors = Dict(
    :Φ=>"Potential",
    :Φ_ext=>"ExtPotential",
    :index=>"ParticleIDs",
    :pos=>"Coordinates",
    :vel=>"Velocities",
    :acc=>"Acceleration",
   )

const snapcols = (:pos, :vel, :acc, :Φ, :Φ_ext, :h, :filename, :header, :m, :index)
const partcols = (:pos, :vel, :acc, :Φ, :Φ_ext, :m, :index)



Base.@kwdef struct Particle
    pos::Vector{F}
    vel::Vector{F}
    acc::Vector{F}
    m::F
    Φ::F = nothing
    Φ_ext::F = nothing
    h::F = NaN
    index::Int
end

function Base.:(==)(p::Particle, q::Particle)
    for sym in partcols
        x = getproperty(p, sym) 
        y = getproperty(q, sym) 
        if any(x .=== y)
            return false
        end
    end
    return true
end


Base.@kwdef struct Snapshot <: AbstractArray{Particle, 1}
    pos::Matrix{F}
    vel::Matrix{F}
    m::F
    acc::Matrix{F} = zeros(size(pos))
    Φ::Vector{F} = zeros(size(pos, 2))
    Φ_ext::Vector{F} = zeros(size(pos, 2))
    header::Dict{String, Any} = make_default_header(size(pos, 2), m)
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
        kwargs[:m] = header["MassTable"][2]
        kwargs[:filename] = filename
    end
    return Snapshot(; kwargs...)
end



# snapshot methods

Base.size(snap::Snapshot) = (length(snap.index),)
Base.IndexStyle(::Type{<:Snapshot}) = IndexLinear()


function iloc(snap::Snapshot, i::Int)
    # return sortperm(snap.index)[i]
end

function Base.getindex(snap::Snapshot, i::Int)
    kwargs = Dict{Symbol, Any}()
    kwargs[:m] = snap.m
    kwargs[:h] = snap.h

    for sym in [:Φ, :Φ_ext, :index]
        kwargs[sym] = getproperty(snap, sym)[i]
    end
    for sym in [:pos, :vel, :acc]
        kwargs[sym] = getproperty(snap, sym)[:, i]
    end
    return Particle(; kwargs...)
end


function Base.getindex(snap::Snapshot, idx::Union{UnitRange, Vector, Colon, BitVector})
    kwargs = Dict{Symbol, Any}()
    kwargs[:m] = snap.m
    kwargs[:h] = snap.h
    kwargs[:header] = snap.header
    kwargs[:filename] = snap.filename

    for sym in [:Φ, :Φ_ext, :index]
        kwargs[sym] = getproperty(snap, sym)[idx]
    end
    for sym in [:pos, :vel, :acc]
        kwargs[sym] = getproperty(snap, sym)[:, idx]
    end
    return Snapshot(; kwargs...)
end


function Base.setindex!(snap::Snapshot, p::Particle, i::Int)
    kwargs = Dict{Symbol, Any}()

    for sym in [:pos, :vel, :Φ, :Φ_ext, :acc]
        getproperty(snap, sym)[i] = getproperty(particle, sym)
    end
    return snap[idx]
end





function Base.copy(snap::Snapshot)
    kwargs = Dict{Symbol, Any}()
    for sym in snapcols
        kwargs[sym] = deepcopy(getproperty(snap, sym))
    end

    return Snapshot(; kwargs...)
end
# high level methods

function write(snap::Snapshot)
    write!(snap.filename, snap)
end

function write!(filename::String, snap::Snapshot)
    h5open(filename, "w") do h5f
        for (var, col) in h5vectors
            val = getproperty(snap, var) # matrix is unneccesary,
            set_vector!(h5f, col, val)
        end

        set_header!(h5f, snap.header)
    end
end






function Base.show(io::IO, p::Particle)
    @printf io "particle at (%4.2f, %4.2f, %4.2f)" p.pos...
    return io
end

function Base.show(io::IO, snap::Snapshot)
    print(io, "<snapshot with $(length(snap)) particles>")
    return io
end


r(p::Particle) = norm(p.pos)
v(p::Particle) = norm(p.vel)


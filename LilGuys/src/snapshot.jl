module SnapshotUtils
export Particle, Snapshot, write!
export r, v

using ..Units
using ..Coordinates
using ..HDF5Utils

using HDF5
using LinearAlgebra: norm
using Printf


const h5scalars = Dict(
    :Φ=>"Potential",
    :Φ_ext=>"ExtPotential",
    :index=>"ParticleIDs"
   )

const h5vectors = Dict(
    :pos=>"Coordinates",
    :vel=>"Velocities",
    :acc=>"Acceleration",
   )

const snapcols = (:pos, :vel, :acc, :Φ, :Φ_ext, :h, :filename, :header, :m, :index)


Base.@kwdef struct Particle
    pos::Point
    vel::Point
    acc::Point
    m::F
    Φ::F = NaN
    Φ_ext::F = NaN
    h::F = NaN
end


Base.@kwdef mutable struct Snapshot <: AbstractArray{Particle, 1}
    pos::Vector{Point}
    vel::Vector{Point}
    acc::Vector{Point}
    Φ::Vector{F}
    Φ_ext::Vector{F}
    index::Vector{Int}
    h::F = NaN
    filename::String
    header::Dict{String, Any}
    m::F
end


function Snapshot(filename::String)
    kwargs = Dict{Symbol, Any}()

    h5open(filename, "r") do h5f
        index = get_vector(h5f, "ParticleIDs")
        perm = sortperm(index)

        for (var, header) in h5scalars
            kwargs[var] = get_vector(h5f, header)[perm]
        end
        for (var, header) in h5vectors
            kwargs[var] = get_vector(h5f, header)[:, perm]
        end

        header = get_header(h5f)
        kwargs[:header] = header
        kwargs[:m] = header["MassTable"][1]
        kwargs[:filename] = filename

    end
    return Snapshot(; kwargs...)
end

# snapshot methods

Base.size(snap::Snapshot) = (length(snap.index),)
Base.IndexStyle(::Type{<:Snapshot}) = IndexLinear()

function Base.getindex(snap::Snapshot, i::Int)
    kwargs = Dict{Symbol, Any}()
    kwargs[:m] = snap.m
    kwargs[:h] = snap.h

    for sym in [:pos, :vel, :Φ, :Φ_ext, :acc]
        kwargs[sym] = getproperty(snap, sym)[i]
    end
    return Particle(; kwargs...)
end


function Base.setindex!(snap::Snapshot, p::Particle, i::Int)
    kwargs = Dict{Symbol, Any}()

    for sym in [:pos, :vel, :Φ, :Φ_ext, :acc]
        getproperty(snap, sym)[i] = getproperty(particle, sym)
    end
    return snap[i]
end


function Base.show(io::IO, snap::Snapshot)
    print(io, "<snapshot with $(length(snap)) particles>")
    return io
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
    write(snap.filename, snap)
end

function write!(filename::String, snap::Snapshot)
    h5open(filename, "w") do h5f
        for (var, col) in h5scalars
            val = getproperty(snap, var)
            set_vector!(h5f, col, val)
        end
        for (var, col) in h5vectors
            val = Matrix(getproperty(snap, var))
            set_vector!(h5f, col, val)
        end

        set_header!(h5f, snap.header)
    end
end



function index_of(snap::Snapshot, id::Int)
    findfirst(x->x==i, snap.index)
end



function Base.show(io::IO, p::Particle)
    @printf io "particle at (%4.2f, %4.2f, %4.2f)" p.pos...
    return io
end


r(p::Particle) = norm(p.pos)
v(p::Particle) = norm(p.vel)

end # module

module SnapshotUtils

using ..Units
using ..Coordinates
using ..HDF5Utils

using HDF5
using Printf



Base.@kwdef struct Particle
    pos::Point
    vel::Point
    m::F
    Φ::F = NaN
    Φ_ext::F = NaN
    h::F = NaN
end


Base.@kwdef mutable struct Snapshot
    pos::Vector{Point}
    vel::Vector{Point}
    Φ::Vector{F}
    Φ_ext::Vector{F}
    index::Vector{Int}
    m::F
    h::F = NaN
end


function Snapshot(filename::String)
    h5f = h5open(filename, "r")
    pos = get_vector(h5f, "Coordinates")
    vel = get_vector(h5f, "Velocities")
    Φ = get_vector(h5f, "Potential")
    Φ_ext = get_vector(h5f, "ExtPotential")
    index = get_vector(h5f, "ParticleIDs")
    
    header = get_header(h5f)
    m = header["MassTable"][1]
    return Snapshot(pos=pos, vel=vel, Φ=Φ, Φ_ext=Φ_ext, index=index, m=m)
end


function Base.getindex(snap::Snapshot, i)
    idx = findfirst(x->x==i, snap.index)
    return Particle(
        pos=snap.pos[idx],
        vel=snap.vel[idx],
        m=snap.m,
        h=snap.h,
        Φ=snap.Φ[idx],
        Φ_ext = snap.Φ_ext[idx]
       )
end


function Base.show(io::IO, snap::Snapshot)
    print(io, "<snapshot with $(length(snap)) particles>")
    return io
end

function Base.show(io::IO, p::Particle)
    @printf io "particle at (%4.2f, %4.2f, %4.2f)" p.pos...
    return io
end

function Base.length(snap::Snapshot)
    return length(snap.index)
end


end

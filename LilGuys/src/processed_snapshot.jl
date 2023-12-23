import Base: @kwdef
using HDF5
using LinearAlgebra: norm
using Printf



"""
A processed snapshot, used during centre-finding and adding stars.
"""
Base.@kwdef mutable struct PSnapshot <: AbstractSnapshot
    snap::Snapshot
    δr::OptVector = nothing
    δv::OptVector = nothing
    w::OptVector = nothing
    cen::FuzzyPhase = nothing
end


function PSnapshot(snap::Snapshot)
    kwargs = Dict{Symbol, Any}()

    for field in fieldnames(Snapshot)
        if getproperty(snap, field) === nothing
            continue
        end
        kwargs[field] = getproperty(snap, field)
    end

    return SSnapshot(;kwargs...)
end

function Base.getproperty(snap::PSnapshot, field::Symbol)
    if field ∈ fieldnames(PSnapshot)
        return getproperty(snap, field)
    else
        return getproperty(snap.snap, field)
    end
end


function Base.copy(snap::PSnapshot)
    kwargs = Dict{Symbol, Any}()
    for field in fieldnames(SSnapshot)
        kwargs[field] = getproperty(snap, field)
    end
    return SSnapshot(;kwargs...)
end


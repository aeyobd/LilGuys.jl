import Base: @kwdef
using HDF5
using LinearAlgebra: norm
using Printf



"""
A snapshot including additionally weights, uncertanties, and an offset.
"""
Base.@kwdef mutable struct FSnapshot <: AbstractSnapshot
    snap::Snapshot
    cen::FuzzyPhase = nothing
    δr::OptVector = nothing
    δv::OptVector = nothing
    w::OptVector = nothing
end

Base.size(snap::FSnapshot) = size(snap.snap)
Base.IndexStyle(::Type{<:FSnapshot}) = IndexLinear()


function FSnapshot(snap::Snapshot)
    kwargs = Dict{Symbol, Any}()
    kwargs[:snap] = copy(snap)
    kwargs[:cen] = FuzzyPhase(zeros(3), zeros(3), 1e5, 1e3)

    return FSnapshot(;kwargs...)
end


function Base.getproperty(snap::FSnapshot, field::Symbol)
    if field ∈ fieldnames(FSnapshot)
        return getfield(snap, field)
    else
        return getfield(snap.snap, field)
    end
end


function Base.copy(snap::FSnapshot)
    kwargs = Dict{Symbol, Any}()
    for field in fieldnames(FSnapshot)
        val = getproperty(snap, field)
        if val === nothing
            ;
        else 
            kwargs[field] = copy(val)
        end
    end
    return FSnapshot(;kwargs...)
end


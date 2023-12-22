using HDF5
using LinearAlgebra: norm
using Printf


Base.@kwdef mutable struct SSnapshot <: AbstractSnapshot
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
    δr::Vector{F} = []
    δv::Vector{F} = []
    w::Vector{F} = []
end

function SSnapshot(snap::Snapshot)
    kwargs = Dict(
        :pos => snap.pos,
        :vel => snap.vel,
        :m => snap.m,
        :acc => snap.acc,
        :Φ => snap.Φ,
        :Φ_ext => snap.Φ_ext,
        :header => snap.header,
        :index => snap.index,
        :filename => snap.filename,
        :h => snap.h,
    )
    return SSnapshot(;kwargs...)
end


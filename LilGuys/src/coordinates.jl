using Printf
import Base: +, *


abstract type CoordinateFrame end

struct ConstVector <: AbstractArray{F, 1}
    value::F
    size::Int
end

function (*)(a::ConstVector, s::F)
    return ConstVector(a.value * s, a.size)
end

function (*)(s::F, a::ConstVector)
    return a * s
end

Base.getindex(v::ConstVector, i::Int) = v.value
Base.show(io::IO, v::ConstVector) = print(io, v.value)
Base.size(v::ConstVector) = (v.size,)
Base.IndexStyle(::Type{<:ConstVector}) = IndexLinear()


@kwdef struct PhasePoint{T} <: CoordinateFrame
    x::F
    y::F
    z::F
    v_x::F = NaN
    v_y::F = NaN
    v_z::F = NaN
end


@kwdef struct SkyCoord{T} <: CoordinateFrame
    ra::F
    dec::F
    distance::F = NaN
    pmra::F = NaN
    pmdec::F = NaN
    radial_velocity::F = NaN
end

const Galactocentric = PhasePoint{:Galactocentric}
const SimPoint = PhasePoint{:SimPoint}
const ICRS_Cartesian = PhasePoint{:ICRS}
const HelioRest_Cartesian = PhasePoint{:HelioRest}

const HelioRest = SkyCoord{:HelioRest}
const ICRS = SkyCoord{:ICRS}
    

function PhasePoint{T}(pos::Vector{F}, vel::Vector{F} = [NaN, NaN, NaN]) where T
    if length(pos) != 3
        error("position must be a 3-vector")
    end
    if length(vel) != 3
        error("velocity must be a 3-vector")
    end
    return PhasePoint{T}(pos..., vel...)
end


function Base.getproperty(pp::PhasePoint, sym::Symbol)
    if sym == :position
        return [pp.x, pp.y, pp.z]
    elseif sym == :velocity
        return [pp.v_x, pp.v_y, pp.v_z]
    elseif sym ∈ fieldnames(PhasePoint)
        return getfield(pp, sym)
    else
        error("$(sym) is not a valid field")
    end
end




function Base.show(io::IO, pp::PhasePoint{T}) where T
    x = pp.position 
    v = pp.velocity
    print(io, "$T point at ")
    @printf io "(%4.2f, %4.2f, %4.2f) kpc, " x...
    @printf io "(%4.2f, %4.2f, %4.2f) km/s" v...
    return io
end


function Base.show(io::IO, coord::SkyCoord{T}) where T
    print(io, "$T at ")
    @printf io "(%4.2f, %4.2f) deg, " coord.ra coord.dec
    @printf io "d = %4.2f kpc, " coord.distance
    @printf io "μ = (%4.2f, %4.2f) mas/yr, " coord.pmra coord.pmdec

    return io
end

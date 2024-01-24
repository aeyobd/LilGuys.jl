using Printf
import Base: +, *


Base.@kwdef struct Particle
    mass::F
    position::Vector{F}
    velocity::Vector{F}
    index::Int
    acceleration::OptVector = nothing
    Φ::F = NaN
    Φ_ext::F = NaN
    h::F = NaN
    weight::F = NaN
end


function approx(x, y; tol=1e-6)
    if x === nothing || y === nothing || x === NaN || y === NaN
        if x !== y
            return false
        end
        return true
    end

    if x == 0 || y == 0
        rerr = @. abs(x - y)
    else
        rerr = @. abs(x - y) / max(abs(x), abs(y))
    end

    if any(rerr .> tol)
        println("x = $x")
        println("y = $y")
        return false
    end
    return true
end


function Base.:(==)(p::Particle, q::Particle)
    for sym in fieldnames(Particle)
        x = getproperty(p, sym) 
        y = getproperty(q, sym) 
        if !approx(x, y)
            return false
        end
    end
    return true
end


function Base.show(io::IO, p::Particle)
    @printf io "particle at (%4.2f, %4.2f, %4.2f)" p.pos...
    return io
end


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

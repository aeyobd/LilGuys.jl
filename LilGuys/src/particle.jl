using Printf
import Base: +, *

OptVector = Union{Vector{F}, Nothing}
OptMatrix = Union{Matrix{F}, Nothing}
OptF = Union{F, Nothing}


Base.@kwdef struct Particle
    m::F
    pos::Vector{F}
    vel::Vector{F}
    index::Int
    acc::OptVector = nothing
    Φ::OptF = nothing
    Φ_ext::OptF = nothing
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

Base.@kwdef struct PhasePoint
    pos::Vector{F}
    vel::Vector{F}
end

Base.@kwdef struct FuzzyPhase
    pos::Vector{F}
    vel::Vector{F}
    δx::F
    δv::F
end

function (+)(p::FuzzyPhase, q::FuzzyPhase)
    return FuzzyPhase(p.pos + q.pos, p.vel + q.vel, p.δx + q.δx, p.δv + q.δv)
end

function Base.copy(p::FuzzyPhase)
    return FuzzyPhase(copy(p.pos), copy(p.vel), p.δx, p.δv)
end

function (*)(a::F, p::FuzzyPhase)
    return FuzzyPhase(a * p.pos, a * p.vel, a * p.δx, a * p.δv)
end


struct ConstVector <: AbstractArray{F, 1}
    val::F
    size::Int
end

Base.getindex(v::ConstVector, i::Int) = v.val
Base.show(io::IO, v::ConstVector) = print(io, v.val)
Base.size(v::ConstVector) = (v.size,)
Base.IndexStyle(::Type{<:ConstVector}) = IndexLinear()

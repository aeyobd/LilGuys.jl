
import Base: +, -, *, copy, getindex

struct Point <: AbstractArray{F, 1}
    x::F
    y::F
    z::F
end

Base.size(p::Point) = (3,)
Base.IndexStyle(::Type{<:Point}) = IndexLinear()

function Base.getindex(p::Point, i::Int)
    if i == 1
        return p.x
    elseif i==2
        return p.y
    elseif i==3
        return p.z
    else
        error(BoundsError("max index is 3, got ", i))
    end
end

function Base.convert(::Type{Point}, p::Vector)
    return Point(p...)
end


Base.length(p::Point) = 3

function (+)(p::Point, q::Point)
    return Point(p.x+q.x, p.y+q.y, p.z+q.z)
end

function (-)(p::Point, q::Point)
    return Point(p.x-q.x, p.y-q.y, p.z-q.z)
end

function (*)(p::Point, a::F)
    return Point(p.x * a, p.y * a, p.z * a)
end

(*)(a::F, p::Point) = p*a

function Base.Vector{Point}(m::Matrix)
    s = size(m)
    if s[1] == 3
        N = s[2]
        return [Point(m[:,i]) for i in 1:N]
    elseif s[1] == 3
        N = s[2]
        return [Point(m[i,:]) for i in 1:N]
    end
end


function Base.Vector(p::Point)
    return [x for x in p]
end
function Base.Matrix(q::Vector{Point})
    return hcat(Vector.(q)...)
end


function Point(v::Vector)
    if length(v) == 3
        return Point(v...)
    end
end

Base.@kwdef struct PhasePoint
    pos::Vector{F}
    vel::Vector{F}
end

Base.@kwdef struct FuzzyPhase
    pos::Point
    vel::Point
    δx::F
    δv::F
end

function (+)(p::FuzzyPhase, q::FuzzyPhase)
    return FuzzyPhase(p.pos + q.pos, p.vel + q.vel, p.δx + q.δx, p.δv + q.δv)
end

function copy(p::FuzzyPhase)
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

using Printf

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


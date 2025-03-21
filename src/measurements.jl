import Base: +, *, ^, -, /, log10, <

"""
    Measurement(val, em, ep)

A type representing an asymmetric error bar with lower and upper errors.
This type propogates errors linearly (worst case).
"""
Base.@kwdef struct Measurement{T<:Real} <: Real
    value::T
    em::T
    ep::T
    note::String = ""
end

function Measurement(value)
    return Measurement(value, zero(value), zero(value))
end

function Measurement{T}(value) where T <: Real
    return Measurement{T}(value, zero(value), zero(value))
end


function Measurement(value, uncertainty, note::String="")
    return Measurement(value, uncertainty, uncertainty, note)
end

function Measurement{T}(value, uncertainty, note::String="") where T<:Real
    return Measurement{T}(value, uncertainty, uncertainty, note)
end

function Measurement(value, em::Real, ep::Real, note::String="")
    return Measurement(value, em, ep, note)
end

function Measurement{T}(value, em::Real, ep::Real) where T <: Real
    return Measurement{T}(value, em, ep, "")
end


function value(x::Measurement)
    return x.value
end

"""
    ci_of(x::Measurement)

Returns the confidence interval of x.
"""
function ci_of(x::Measurement)
    return x.em, x.ep
end

function (+)(a::Measurement, b::Measurement)
    return Measurement(a.value+b.value, a.em+b.em, a.ep+b.ep, a.note * "; " * b.note)
end

function (-)(a::Measurement, b::Measurement)
    return Measurement(a.value-b.value, a.em+b.em, a.ep+b.ep, a.note * "; " * b.note)
end


function log10(a::Measurement)
    m = log10(a.value)
    l = log10(a.value - a.em)
    h = log10(a.value + a.em)
    return Measurement(m, m-l, h-m, a.note)
end


function (^)(b::Real, a::Measurement)
    m = b^(a.value)
    l = b^(a.value - a.em)
    h = b^(a.value + a.em)
    return Measurement(m, m-l, h-m, a.note)
end


function (+)(x::Measurement, y::Real)
    return Measurement(x.value+y, x.em, x.ep, x.note)
end

function (-)(x::Measurement, y::Real)
    return Measurement(x.value-y, x.em, x.ep, x.note)
end

function (*)(x::Measurement, y::Real)
    return Measurement(x.value*y, x.em*y, x.ep*y, x.note)
end

function (*)(x::Measurement, y::Measurement)
    m = x.value*y.value
    l = (x.value - x.em) * (y.value - y.em)
    u = (x.value + x.em) * (y.value + y.em)

    return Measurement(m, m-l, u-m)
end

function (/)(x::Measurement, y::Real)
    return Measurement(x.value/y, x.em/y, x.ep/y, x.note)
end

function (/)(x::Measurement, y::Measurement)

    m = x.value/y.value
    l = (x.value - x.em) / (y.value + y.em)
    u = (x.value + x.em) / (y.value - y.em)

    return Measurement(m, m-l, u-m)
end

Base.promote_rule(::Type{<:Real}, ::Type{Measurement{T}}) where {T<:Real} = Measurement{T}


function Base.print(io::IO, x::Measurement)
    print(io, "$(x.value) + $(x.ep) - $(x.em)")
end

function Base.isapprox(a::Measurement, b::Measurement; kwargs...)
    return (
        isapprox(a.value, b.value; kwargs...) && 
        isapprox(a.em, b.em; kwargs...) && 
        isapprox(a.ep, b.ep; kwargs...)
       )
end

function Base.isfinite(a::Measurement)
    return (isfinite.(a.value) && isfinite.(a.em) && isfinite.(a.ep))
end

function (<)(a::Measurement, b::Measurement)
    return a.value < b.value
end
function (<)(a::Measurement, b::Real)
    return a.value < b
end

function (-)(a::Measurement)
    return Measurement(-a.value, a.em, a.ep, a.note)
end

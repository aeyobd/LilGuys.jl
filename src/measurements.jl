import Base: +, *, ^, -, /, log10, <

"""
    Measurement(middle, lower, upper, kind)

A type representing a measurement with lower and upper errors.
This type propogates errors linearly (worst case).

`kind` is a note to self. Recommended to be formatted like

- normal
- 0.16_quantile or 0.05_quantile
- 0.16_mode_hdi

Will warn if annotations do not match when combining measurements.
"""
Base.@kwdef struct Measurement{T<:Real} <: Real
    middle::T
    lower::T
    upper::T
    kind::String = ""
end


function Measurement(middle)
    return Measurement(middle, zero(middle), zero(middle))
end


function Measurement{T}(middle) where T <: Real
    return Measurement{T}(middle, zero(middle), zero(middle))
end


function Measurement(middle, uncertainty, kind::String="")
    return Measurement(middle, uncertainty, uncertainty, kind)
end


function Measurement{T}(middle, uncertainty, kind::String="") where T<:Real
    return Measurement{T}(middle, uncertainty, uncertainty, kind)
end


function Measurement(middle, lower::Real, upper::Real)
    return Measurement(middle, lower, upper, "")
end


function Measurement{T}(middle, lower::Real, upper::Real) where T <: Real
    return Measurement{T}(middle, lower, upper, "")
end


"""
    middle(x::Measurement)

Return the middle/central value of x.
"""
function middle(x::Measurement)
    return x.middle
end


"""
    error_interval(x::Measurement)

Returns the error interval on x, i.e. the lower and upper errors
as a tuple of positive numbers.
"""
function error_interval(x::Measurement)
    return x.lower, x.upper
end


"""
    lower_error(x::Measurement)

Return the lower error/uncertainty/credible range deviation of x.
"""
function lower_error(x::Measurement)
    return x.lower
end


"""
    upper_error(x::Measurement)

Return the upper error/uncertainty/credible range deviation of x.
"""
function upper_error(x::Measurement)
    return x.upper
end


function common_kind(a::Measurement, b::Measurement)
    if a.kind == b.kind
        kind = a.kind
    else
        @warn "measurement kinds do not match" 
        kind = ""
    end

    return kind

end


function (+)(a::Measurement, b::Measurement)
    return Measurement(a.middle+b.middle, a.lower+b.lower, a.upper+b.upper, common_kind(a, b))
end


function (-)(a::Measurement, b::Measurement)
    return Measurement(a.middle-b.middle, a.lower+b.lower, a.upper+b.upper, common_kind(a, b))
end


function log10(a::Measurement)
    m = log10(a.middle)
    l = log10(a.middle - a.lower)
    h = log10(a.middle + a.lower)
    return Measurement(m, m-l, h-m, a.kind)
end


function (^)(b::Real, a::Measurement)
    m = b^(a.middle)
    l = b^(a.middle - a.lower)
    h = b^(a.middle + a.lower)
    return Measurement(m, m-l, h-m, a.kind)
end


function (+)(x::Measurement, y::Real)
    return Measurement(x.middle+y, x.lower, x.upper, x.kind)
end


function (-)(x::Measurement, y::Real)
    return Measurement(x.middle-y, x.lower, x.upper, x.kind)
end


function (*)(x::Measurement, y::Real)
    return Measurement(x.middle*y, x.lower*y, x.upper*y, x.kind)
end


function (*)(x::Measurement, y::Measurement)
    m = x.middle*y.middle
    l = (x.middle - x.lower) * (y.middle - y.lower)
    u = (x.middle + x.lower) * (y.middle + y.lower)

    return Measurement(m, m-l, u-m, common_kind(x, y))
end


function (/)(x::Measurement, y::Real)
    return Measurement(x.middle/y, x.lower/y, x.upper/y, x.kind)
end


function (/)(x::Measurement, y::Measurement)
    m = x.middle/y.middle
    l = (x.middle - x.lower) / (y.middle + y.lower)
    u = (x.middle + x.lower) / (y.middle - y.lower)

    return Measurement(m, m-l, u-m, common_kind(x, y))
end


function Base.promote_rule(::Type{<:Real}, ::Type{Measurement{T}}) where {T<:Real}
    return Measurement{T}
end


function Base.print(io::IO, x::Measurement)
    print(io, "$(x.middle) + $(x.upper) - $(x.lower)")
end


function Base.isapprox(a::Measurement, b::Measurement; kwargs...)
    return (
        isapprox(a.middle, b.middle; kwargs...) && 
        isapprox(a.lower, b.lower; kwargs...) && 
        isapprox(a.upper, b.upper; kwargs...)
       )
end


function Base.isfinite(a::Measurement)
    return (isfinite.(a.middle) && isfinite.(a.lower) && isfinite.(a.upper))
end


function (<)(a::Measurement, b::Measurement)
    return a.middle < b.middle
end


function (<)(a::Measurement, b::Real)
    return a.middle < b
end


function (-)(a::Measurement)
    return Measurement(-a.middle, a.lower, a.upper, a.kind)
end

import Base: *


"""
    ConstVector(value, size)

Memory efficient vector which is constant (useful for type-enforced vectors, e.g. unit weights...)
"""
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


function Base.getindex(v::ConstVector, i::Int)
    if i < 1 || i > v.size
        throw(BoundsError(v, i))
    end
    return v.value
end

Base.show(io::IO, v::ConstVector) = print(io, v.value)
Base.size(v::ConstVector) = (v.size,)
Base.IndexStyle(::Type{<:ConstVector}) = IndexLinear()


"""
general method to convert a struct to a dictionary
"""
function struct_to_dict(S, split_errors=true)
    d = Dict(key=>getfield(S, key) for key in fieldnames(typeof(S)))

    if split_errors
        for (key, val) in d
            if val isa AbstractArray{<:Measurement}
                m = LilGuys.middle.(val)
                em = LilGuys.lower_error.(val)
                ep = LilGuys.upper_error.(val)
                d[key] = m
                d[Symbol(string(key) * "_em")] = em
                d[Symbol(string(key) * "_ep")] = ep
            end
        end
    end

    return d
end

function collapse_errors(d::Dict)
    d_new = deepcopy(d)
    ks = keys(d) 
    for (key, val) in d
        if [key * "_em", key*"_ep"] ⊆ ks
            em = pop!(d_new, key*"_em")
            ep = pop!(d_new, key*"_ep")
            d_new[key] = Measurement.(d_new[key], em, ep)
        end

    end

    return d_new
end


"""
general method to convert a dictionary to a named tuple
"""
function dict_to_tuple(D)
    return NamedTuple((Symbol(key), value) for (key, value) in D)
end



"""
    randu(low, high[, size...])

Returns a random number between low and high.
"""
function randu(low::Real, high::Real, args...)
    return low .+ (high - low) * rand(args...)
end


"""
    rand_unit(N::Integer)

Returns a 3xN matrix of random unit vectors.
"""
function rand_unit(N::Integer=1)
    # generate a random vector
    x = randn(3, N)
    x_norm = reshape(radii(x), 1, N)

    return x ./ x_norm
end




"""
    gradient(y[, x])

computes the gradient (dy/dx) of a 2D function at the point (x, y).
assumes that x are sorted.
Returns a vector same length of x with endpoints using linear approximation.
Uses the 2nd order central difference method alla numpy.gradient.
"""
function gradient(y::AbstractVector{T}, x::AbstractVector) where T<:Real
    x = x
    y = y
    N = length(x)

    if N < 2
        return fill(NaN, N)
    end

    grad = Vector{T}(undef, N)

    grad[1] = (y[2] - y[1]) / (x[2] - x[1])
    grad[end] = (y[end] - y[end-1]) / (x[end] - x[end-1])
    for i in 2:(N-1)
            hs = x[i] - x[i-1]
            hd = x[i+1] - x[i]

            numerator = hs^2 * y[i+1] + (hd^2 - hs^2) * y[i] - hd^2*y[i-1]
            denom = hd*hs*(hd + hs)
            grad[i] = numerator/denom
    end
    return grad
end



"""
    gradient(y)
computes the gradient
"""
function gradient(y::AbstractVector{T}) where T<:Real
    x = collect(1.0:length(y))
    return gradient(y, x)
end



"""
Represents a linear interpolator. Construct via `lerp` for now,
expects `x` to be a sorted list to effieintly returl interpolated values.
"""
struct LinInterp{T<:Real}
    x::Vector{T}
    y::Vector{T}
end


function (l::LinInterp)(x::Real)
    xs = l.x
    ys = l.y

    if !isfinite(x)
        return NaN
    end
    if x < xs[1]
        return ys[1]
    elseif x > xs[end]
        return ys[end]
    end
    i = searchsortedfirst(xs, x)
    if i <= 1
        return ys[1]
    elseif i > length(xs)
        return ys[end] # TODO: test
    elseif 1 < i <= length(xs)
        x1, x2 = xs[i-1], xs[i]
        y1, y2 = ys[i-1], ys[i]
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    else
        return NaN # TODO: test
    end
end


"""
Returns a linear interpolation of the given xs and ys.
Is truncated and will return the first or last value if x is outside the range.
"""
function lerp(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    if length(xs) != length(ys)
        throw(ArgumentError("xs and ys must have the same length, got $(length(xs)) and $(length(ys))")) 
    end
    ys = collect(copy(ys))
    xs = collect(copy(xs))
    ys = ys[sortperm(xs)]
    xs = sort(xs)

    return LinInterp(xs, ys)
end



function normal_cdf(x::Real, μ::Real, σ::Real)
    z = (x - μ) / σ
    return normal_cdf(z)
end


function normal_cdf(z::Real)
    return 1/2 * (1 + erf(z / √2) )
end

function gaussian(z::Real)
    return exp(-z^2 / 2) / √(2π)
end

function gaussian(x::Real, μ::Real, σ::Real)
    z = (x - μ) / σ
    return gaussian(z) / σ
end


"""
    logistic(z)

Computes the logistic function at the point z.
```math
f(z) = \\frac{1}{1 + e^{-z}} 
```
"""
function logistic(z::Real)
    return 1 / (1 + exp(-z))
end


"""
    logit(p)

Computes the logit function at the point p.
```math
f(p) = \\log\\left(\\frac{p}{1-p}\\right)
````
"""
function logit(p::Real) # TODO: test
    return log(p / (1-p))
end



"""
    centroid(x_vec[, weights])

Computes the centroid of a 3xN matrix of points, returning the mean and standard deviation.
"""
function centroid(x_vec::AbstractMatrix{T}) where T<:Real
    cen = mean(x_vec, dims=2)
    return cen[:, 1]
end



function centroid(x::AbstractMatrix{T1}, weights::AbstractVector{T2}) where T1 <: Real where T2<:Real
    N = size(x, 2)
    w = reshape(weights, :, 1) ./ sum(weights)
    cen = (x * w)
    return cen[:, 1]
end


"""
    centroid_err(x_vec[, weights])

Computes the standard error of the centroid of a 3xN matrix of points.
"""
function centroid_err(x::Matrix{T}) where T<:Real
    N = size(x, 2)
    σs = std(x, dims=2) / sqrt(N) # standard error of mean
    return sqrt(mean(σs.^2))
end



function centroid_err(x::AbstractMatrix{T}, weights::AbstractVector) where T<:Real
    N = size(x, 2)
    if N <= 1
        return NaN
    end
    c = centroid(x, weights)
    w = reshape(weights, :, 1) ./ sum(weights)
    s = mean((x .- c).^2 * w)
    err = sqrt(s) / sqrt(N-1)
    return err
end





"""
    sample_surface_density(f::Function, N::Integer = 1; log_R=nothing)

Randomly draws N samples from a Density profile given by the finction f.
"""
function sample_surface_density(f::Function, N::Integer = 1; log_R=nothing)
    if log_R == nothing
        log_R = LinRange(-5, 5, 1000)
    end

    R = exp10.(log_R)
    Σ = f.(R)

    if any(Σ .< 0)
        throw(ArgumentError("Σ must be positive"))
    end
    if any(.! isfinite.(Σ))
        throw(ArgumentError("Σ must be finite"))
    end


    M = cumsum(Σ .* π .* R .* gradient(R))
    M = M ./ M[end]

    l = lerp([0; M], [0; R])

    probs = rand(N)
    return l.(probs)
end


"""
    sample_density(f::Function, N::Integer = 1; log_r=nothing)

Randomly draws N samples from a Density profile given by the finction f.
"""
function sample_density(f::Function, N::Integer = 1; log_r=nothing)
    if log_r == nothing
        log_r = LinRange(-5, 5, 10000)
    end

    r = 10 .^ log_r
    ρ = f.(r)

    if any(ρ .< 0)
        throw(ArgumentError("ρ must be positive"))
    end
    if any(.! isfinite.(ρ))
        throw(ArgumentError("ρ must be finite"))
    end


    dr = gradient(r)
    M = cumsum(@. ρ * 4π*r^2 * dr)
    M = M ./ M[end]

    l = lerp([0; M], [0; r])

    probs = rand(N)
    return l.(probs)
end


"""
    @assert_same_size x y

Asserts that x and y are the same size
"""
macro assert_same_size(x, y)
    return quote
        if size($(esc(x))) != size($(esc(y)))
            throw(ArgumentError("$(string($(QuoteNode(x)))) and $(string($(QuoteNode(y)))) must have the same size"))
        end
    end
end


macro assert_3vector(x)
    return quote
        local size_x = size($(esc(x)))
        if size_x[1] != 3
            throw(ArgumentError("$(string($(QuoteNode(x)))) must be a 3-vector, got $(size_x)"))
        end
    end
end



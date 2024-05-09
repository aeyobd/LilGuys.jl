import StatsBase: mean, std, percentile
import SpecialFunctions: erf



"""
    randu(low, high[, size...])

Returns a random number between low and high.
"""
function randu(low::Real, high::Real, args...)
    return low .+ (high - low) * rand(args...)
end


"""
    rand_unit(N)

Returns a 3xN matrix of random unit vectors.
"""
function rand_unit(N::Int=1)
    # generate a random vector
    x = randn(3, N)
    x_norm = reshape(calc_r(x), 1, N)

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
    midpoint(x)

Computes the midpoint of each element of a vector (for example, bin centres).
"""
function midpoint(x::AbstractVector{T}) where T<:Real
    return (x[1:end-1] + x[2:end]) / 2
end




"""
Returns a linear interpolation of the given xs and ys
"""
function lerp(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    if length(xs) != length(ys)
        throw(ArgumentError("xs and ys must have the same length, got $(length(xs)) and $(length(ys))"))
    end

    return function(x::Real)
        if x < xs[1]
            return ys[1]
        elseif x > xs[end]
            return ys[end]
        end
        i = searchsortedfirst(xs, x)
        if i <= 1
            return ys[1]
        elseif i > length(xs)
            return ys[end]
        elseif 1 < i <= length(xs)
            x1, x2 = xs[i-1], xs[i]
            y1, y2 = ys[i-1], ys[i]
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        else
            return NaN
        end
    end

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
function logit(p::Real)
    return log(p / (1-p))
end






"""
calculates equal number bins over the array x with n values per bin.
"""
function make_equal_number_bins(x, n)
    return [percentile(x, i) for i in LinRange(0, 100, n+1)]
end


"""
    calc_histogram(x, bins; weights=Nothing)

Computes the histogram of a vector x with respect to the bins with optional weights. Returns the bin edges and the histogram.
"""
function calc_histogram(x::AbstractVector, bins::AbstractVector; weights=Nothing)
    if weights == Nothing
        weights = ones(Int64, length(x))
    end
    Nbins = length(bins) - 1
    hist = zeros(Nbins)
    for i in 1:Nbins
        idx = (x .>= bins[i]) .& (x .< bins[i+1])
        hist[i] = sum(weights[idx])
    end
    return bins, hist
end


function calc_histogram(x::AbstractVector, bins::Int=20; weights=Nothing, xlim=Nothing)
    if xlim == Nothing
        xlim = (minimum(x[isfinite.(x)]), maximum(x[isfinite.(x)]))
    end

    return calc_histogram(x, LinRange(xlim[1], xlim[2], bins); weights=weights)
end




"""
    centroid(x_vec[, weights])

Computes the centroid of a 3xN matrix of points, returning the mean and standard deviation.
"""
function centroid(x_vec::Matrix{T}) where T<:Real
    cen = mean(x_vec, dims=2)
    return cen[:, 1]
end



function centroid(x::Matrix{T}, weights::AbstractVector{T}) where T<:Real
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



function centroid_err(x::Matrix{T}, weights::AbstractVector{T}) where T<:Real
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


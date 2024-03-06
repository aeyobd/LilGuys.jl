import StatsBase: mean, std


"""
    centroid(x_vec[, weights])

Computes the centroid of a 3xN matrix of points, returning the mean and standard deviation.
"""
function centroid(x_vec::Matrix{T}) where T<:Real
    cen = mean(x_vec, dims=2)
    return cen[:, 1]
end



function centroid(x::Matrix{T}, weights::Vector{T}) where T<:Real
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



function centroid_err(x::Matrix{T}, weights::Vector{T}) where T<:Real
    N = size(x, 2)
    c = centroid(x, weights)
    w = reshape(weights, :, 1) ./ sum(weights)
    s = mean((x .- c).^2 * w)
    err = sqrt(s) / sqrt(N-1)
    return err
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



function normal_cdf(x::Real, μ::Real, σ::Real)
    z = (x - μ) / σ
    return normal_cdf(z)
end


function normal_cdf(z::Real)
    return 0.5 * (1 + erf(z / sqrt(2)))
end


"""
    gradient(x, y)

computes the gradient of a 2D function at the point (x, y).
assumes that x are sorted.
"""
function gradient(x::AbstractVector{T}, y::AbstractVector{T}) where T<:Real
	x = x
	y = y
	N = length(x)

	grad = Vector{Float64}(undef, N)

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
    midpoint(x)

Computes the midpoint of each element of a vector (for example, bin centres).
"""
function midpoint(x::AbstractVector{T}) where T<:Real
    return (x[1:end-1] + x[2:end]) / 2
end



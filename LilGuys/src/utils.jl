
"""
The mean of a vector of numbers.
"""
function mean(x::Array{T}) where T<:Real
    N = length(x)
    return sum(x) / N
end

"""
The variance of a vector of numbers.
"""
function var(x::Array{T}) where T<:Real
    N = length(x)
    μ = mean(x)
    return sum((x .- μ).^2) / (N - 1)
end

"""
The standard deviation of a vector of numbers.
"""
function std(x::Array{T}) where T<:Real
    return sqrt(var(x))
end


"""
Computes the centroid of a 3xN matrix of points, returning the mean and standard deviation.
"""
function centroid(x::Matrix{T}) where T<:Real
    c = sum(x, dims=2) / size(x, 2)
    c = c[:, 1]
    variance = mean(sum((x .- c).^2, dims=2) / size(x, 2))
    return c, sqrt(variance)
end

function centroid(x::Matrix{T}, weights::Vector{T}) where T<:Real
    w = reshape(weights, :, 1) ./ sum(weights)
    c = (x * w)[:, 1]
    variance = mean((x .- c).^2 * w)
    err = sqrt(variance)
    return c, err
end


"""
a random unit vector
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


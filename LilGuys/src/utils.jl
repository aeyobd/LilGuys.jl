
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
    N = size(x, 2)
    
    cen = sum(x, dims=2) / N
    return cen[:, 1]
end

function centroid_err(x::Matrix{T}) where T<:Real
    c = centroid(x)
    s = sum((x.-c).^2)
    err = sqrt(s) / √N / sqrt(N-1) # √N for sample error, √N-1 for varience
    return err
end


function centroid(x::Matrix{T}, weights::Vector{T}) where T<:Real
    N = size(x, 2)
    w = reshape(weights, :, 1) ./ sum(weights)

    cen = (x * w)
    return cen[:, 1]
end

function centroid_err(x::Matrix{T}, weights::Vector{T}) where T<:Real
    c = centroid(x, weights)
    s = sum((x .- c).^2 * w)
    err = sqrt(s) / sqrt(N-1)
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


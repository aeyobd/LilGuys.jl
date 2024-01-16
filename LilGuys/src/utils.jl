import StatsBase: mean, std


"""
Computes the centroid of a 3xN matrix of points, returning the mean and standard deviation.
"""
function centroid(x_vec::Matrix{T}) where T<:Real
    cen = mean(x_vec, dims=2)
    return cen[:, 1]
end


function centroid_err(x::Matrix{T}) where T<:Real
    N = size(x, 2)
    σs = std(x, dims=2) / sqrt(N) # standard error of mean
    return sqrt(mean(σs.^2))
end


function centroid(x::Matrix{T}, weights::Vector{T}) where T<:Real
    N = size(x, 2)
    w = reshape(weights, :, 1) ./ sum(weights)
    cen = (x * w)
    return cen[:, 1]
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

function gradient(x::Vector{T}) where T<:Real

end

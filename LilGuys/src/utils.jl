import StatsBase: mean, std

# """
# The mean of a vector of numbers.
# """
# function mean(x::Array{T}) where T<:Real
#     N = length(x)
#     return sum(x) / N
# end
# 
# 
# """
# The variance of a vector of numbers.
# """
# function var(x::Array{T}) where T<:Real
#     N = length(x)
#     μ = mean(x)
#     return sum((x .- μ).^2) / (N - 1)
# end
# 
# 
# """
# The standard deviation of a vector of numbers.
# """
# function std(x::Array{T}) where T<:Real
#     return sqrt(var(x))
# end
# 

"""
Computes the centroid of a 3xN matrix of points, returning the mean and standard deviation.
"""
function centroid(x_vec::Matrix{T}) where T<:Real
    cen = mean(x_vec, dims=2)
    return cen[:, 1]
end

function centroid_err(x::Matrix{T}) where T<:Real
    N = size(x, 2)
    L = length(x)
    return std(x) * sqrt(L - 1) / sqrt(N * (N-1))
end


function centroid(x::Matrix{T}, weights::Vector{T}) where T<:Real
    N = size(x, 2)
    w = reshape(weights, :, 1) ./ sum(weights)
    cen = (x * weights)
    return cen[:, 1]
end

function centroid_err(x::Matrix{T}, weights::Vector{T}) where T<:Real
    N = size(x, 2)
    c = centroid(x, weights)
    s = sum((x .- c).^2 * weights)
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

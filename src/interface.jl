"""
This module provides an interface to external packages.
The only packages I would like in the main library are part of base julia.
Otherwise even simple import and reexport packages from here is ideal, 
especially as Julia is still young.
"""
module Interface

export mean, std, midpoints, weights, StatsBase, quantile, varience
export erf, expinti

export integrate, curve_fit, find_zero
export ±, uncertainty, value, Measurement

export histogram, effective_sample_size, bins_default, bins_equal_number, bins_both, bins_equal_width
export default_bin_width, default_n_per_bin
export DataFrame

# TODO: wrap DataFrame around tables.jl interface and only export `Table` and minimal functions :)
#
# TODO: move io functions here as well 



import StatsBase as sb
import SpecialFunctions: erf, expinti

import QuadGK: quadgk
import LsqFit

using Measurements: ±, value, uncertainty, Measurement
import DensityEstimators

import Roots: find_zero
import DataFrames: DataFrame


midpoints = sb.midpoints

"""
  curve_fit(model, xdata, ydata, p0) -> fit
  curve_fit(model, xdata, ydata, weights, p0) -> fit

Fit the model to the data using the initial guess `p0`. weights are optional. Returns the best fit (Least Squares) parameters and the estimated covariance matrix.
"""
function curve_fit(args...; kwargs...)
    fit = LsqFit.curve_fit(args...; kwargs...)

    if !fit.converged
        @warn "LSq fit did not converge"
    end

    return fit.param, LsqFit.vcov(fit)
end


"""
    integrate(f, a, b; kwargs...)

Integrate the function `f` from `a` to `b`.
A thin wrapper around `quadgk` and prints to @info if the error is > 1e-4.
"""
function integrate(f, a, b; kwargs...)
    val, err = quadgk(f, a, b; kwargs...)
    if err / val > 1e-4
        @info "Warning: relative integration error is large: $(err / val)"
    end

    return val
end


"""
    histogram(x, bins=nothing; weights=nothing, normalization=:none)

Compute the histogram of `x` with `bins` bins, returning the bins, values, and poisson uncertainties.
If `bins` is not provided, it defaults to the rule provided by `default_bins`.
If `weights` are provided, they are used to weight the histogram.
If `normalization` is `:none`, the histogram is normalized to the total number of observations.
"""
function histogram(x, bins=bins_default; 
        weights=nothing, normalization=:none, kwargs...)

    x, weights = filter_nans(x, weights)

    if bins === nothing
        bins = bins_default
    end
    if bins isa Function
        bins = bins(x, weights; kwargs...)
        @info "Using default bins of size = $(length(bins))"
    end

    h = DensityEstimators.histogram(x, bins, weights=weights, normalization=normalization)

    return h.bins, h.values, h.err
end


"""
    filter_nans(x, weights=nothing)

Filter out NaNs from x and weights. If weights are provided, filter out NaNs from
both x and weights. Returns the filtered x and weights.
"""
function filter_nans(x, weights=nothing)
    filt = isfinite.(x)
    if weights !== nothing
        filt .&= isfinite.(weights)
        weights = weights[filt]
    end
    x = x[filt]

    if sum(.!filt) > 0
        @info "Filtered out $(sum(.!filt)) NaNs out of $(length(x))"
    end

    return x, weights
end

"""
    bins_default(x, weights=nothing)

The default setting for bins, returns the equal width bins.
"""
function bins_default(x, weights)
    return bins_equal_width(x, weights)
end


"""
    bins_equal_width(x, weights=nothing)

Returns a vector of bins
"""
function bins_equal_width(x, weights; bin_width=nothing)
    x, weights = filter_nans(x, weights)
    if bin_width === nothing
        bin_width = default_bin_width(x, weights)
    end

    γ = 0.5
    x_l = minimum(x) - γ * bin_width 
    x_u = maximum(x) +  bin_width

    return x_l:bin_width:x_u
end



"""
    bins_equal_number(x, weights; num_bins=nothing)

Returns bins with an equal number of observations per bin
"""
function bins_equal_number(x, weights; num_per_bin=nothing)
    x, weights = filter_nans(x, weights)

    if num_per_bin === nothing
        num_per_bin = default_n_per_bin(x, weights)
        @info "Using $num_per_bin observations per bins"
    end

    num_bins = ceil(Int, length(x) / num_per_bin)
    q = LinRange(0, 1, num_bins + 1)

    bins = quantile(x, q)
    
    bins = unique(bins) 
    @assert length(bins) > 1 "Not enough bins, have $bins"
    return bins
end


"""
    bins_both(x, weights; bin_width=nothing, num_per_bin=nothing)

Computes bins where each bin is of width at least `bin_width` and contains at least `num_per_bin` observations. Returns an array of type x containing the bins.
"""
function bins_both(x, weights; bin_width=nothing, num_per_bin=nothing)
    x, weights = filter_nans(x, weights)

    if bin_width === nothing
        bin_width = default_bin_width(x, weights)
    end
    if num_per_bin === nothing
        num_per_bin = ceil(Int, default_n_per_bin(x, weights))
    end


    bins = empty(x)
    x = sort(x)

    x_i = minimum(x) - 0.5*bin_width
    push!(bins, x_i)
    N = length(x)

    while x_i < maximum(x) - bin_width
        x_binwidth = x_i + bin_width
        idx = searchsortedfirst(x, x_i)

        if idx + num_per_bin >= N # we are out of data
            break
        end
        x_number = x[idx + num_per_bin]
        x_i = max(x_binwidth, x_number)

        push!(bins, x_i)
    end

    if (sum(x .> bins[end]) <= num_per_bin) || (maximum(x) - bins[end] < bin_width)
        bins[end] = maximum(x) + 0.5*bin_width
    else
        push!(bins, maximum(x) + 0.5*bin_width)
    end

    return bins
end


"""
    default_bin_width(x, weights=nothing)

Calculates the Freedman-Diaconis rule for bin size.
"""
function default_bin_width(x, weights=nothing)
    N = effective_sample_size(x, weights)

    iqr = quantile(x, 0.75) - quantile(x, 0.25)
    bin_width = 2 * iqr / N^(1/3)
    return bin_width
end


"""
    default_n_per_bin(x, weights=nothing)

Calculates the recommended number of observations per bin.
"""
function default_n_per_bin(x, weights=nothing)
    N = effective_sample_size(x, nothing) # ignore weight for now

    num_per_bin = 2N^(2/5)
    return num_per_bin * length(x) / N
end


"""
    mean(x; kwargs...)

Returns the arithmatic mean of x. Kwargs may be dims
"""
function mean(x; kwargs...)
    return sb.mean(x; kwargs...)
end

"""
    mean(x, w; kwargs...)

Returns the weighted mean of x with weights w. Kwargs may be dims
"""
function mean(x, w; kwargs...)
    return sb.mean(x, sb.weights(w); kwargs...)
end


"""
    std(x; kwargs...)

Returns the (bias corrected) standard deviation of a vector of real numbers.
"""
function std(x; kwargs...)
    return sb.std(x; corrected=true, kwargs...)
end

"""
    std(x, w; kwargs...)

Returns the weighted standard deviation of x with weights w
"""
function std(x, w; kwargs...)
    return sb.std(x, sb.weights(w); corrected=false, kwargs...)
end

"""
    variance(x; kwargs...)

Returns the (bias corrected) variance of a vector of real numbers.
"""
function variance(x; kwargs...)
    return sb.var(x; corrected=true, kwargs...)
end

"""
    variance(x, w; kwargs...)

Returns the weighted variance of x with weights w
"""
function variance(x, w; kwargs...)
    return sb.var(x, sb.weights(w); corrected=false, kwargs...)
end




"""
    quantile(x, p)

Returns the pth quantile of x
"""
function quantile(x, p)
    return sb.quantile(x, p)
end

"""
    quantile(x, w, p)

Returns the pth quantile of x with weights w
"""
function quantile(x, w, p)
    ws = sb.weights(w)
    return sb.quantile(x, sb.weights(w), p)
end


@doc raw"""
    effective_size(data, weights)

Computes the effective size of a set of weights. If all weights are equal, than
the effective sample size is simply the number of observations (the length of
the weights). However, for more variable weight distributions, the effective
sample size will decrease.

The equation (Kish) is given by
```math
n_{eff} = \frac{ \left( \sum w \right)^2 }{ \sum w^2 } 
```
"""
function effective_sample_size(weights::AbstractVector{<:Real})
    return sum(weights)^2 / sum(weights .^ 2)
end


function effective_sample_size(data::AbstractVector{<:Real}, weights::Nothing)
    return length(data)
end

function effective_sample_size(data::AbstractVector{<:Real}, weights::AbstractVector{<:Real})
    return effective_sample_size(weights)
end

end

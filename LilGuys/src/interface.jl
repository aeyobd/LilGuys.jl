"""
This module provides an interface to external packages.
The only packages I would like in the main library are part of base julia.
Otherwise even simple import and reexport packages from here is ideal, 
especially as Julia is still a baby.
"""
module Interface

export mean, std, percentile, midpoints, weights, StatsBase, quantile
export erf, expinti

export integrate, curve_fit, find_zero
export ±, uncertainty, value, Measurement

export histogram
export DataFrame



import StatsBase: mean, std, percentile, midpoints, StatsBase, weights, quantile
import SpecialFunctions: erf, expinti

import QuadGK: quadgk
import LsqFit: curve_fit

using Measurements: ±, value, uncertainty
import DensityEstimators

import Roots: find_zero
import DataFrames: DataFrame



"""
    integrate(f, a, b; kwargs...)

Integrate the function `f` from `a` to `b`.
A thin wrapper around `quadgk` and prints to @info if the error is > 1e-4.
"""
function integrate(f, a, b; kwargs...)
    val, err = quadgk(f, a, b; kwargs...)
    if err > 1e-4
        @info "Warning: integration error is large: $err"
    end

    return val
end


"""
    histogram(x, bins=nothing; weights=nothing, normalization=:none)

Compute the histogram of `x` with `bins` bins.
If `bins` is not provided, it defaults to 20.
If `weights` are provided, they are used to weight the histogram.
"""
function histogram(x, bins=nothing; weights=nothing, normalization=:none)
    if bins === nothing
        bins = 20
        @info "Using default bins=20"
    end

    h = DensityEstimators.histogram(x, bins, weights=weights, normalization=:none)

    return h.bins, h.values, h.err
end




# TODO: wrap DataFrame around tables.jl interface and only export `Table` and minimal functions :)
#
# TODO: move io functions here as well 

end
#using Glob

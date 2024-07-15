import Arya
import StatsBase: percentile


"""

sorts a snapshot by radius from 0
"""
function sort_by_r(snap::Snapshot)
    return snap[sortperm(calc_r(snap))]
end



"""
    calc_m_hist(r, r_bins[, masses])

Calculates the density profile given a set of particles located at `r` with masses `masses` by binning into `r_bins`.
"""
function calc_ρ_hist(r::AbstractVector{T}, bins::AbstractVector; weights=nothing) where T <: Real
    if weights == nothing
        weights = ones(length(r))
    end

    counts = Arya.histogram(r, bins, weights=weights, normalization=:none).values

    Vs = 4π/3 * diff(bins .^ 3)
    return bins, counts ./ Vs
end


function calc_ρ_hist(r::AbstractVector{T}, bins::Int; weights=nothing, equal_width=false) where T <: Real
    if equal_width
        x1 = minimum(r)
        x2 = maximum(r)
        x = LinRange(x1, x2, bins)
        bins = 10 .^ x
    else
        bins = percentile(r, LinRange(0, 100, bins+1))
    end
    return calc_ρ_hist(r, bins; weights=weights)
end


function calc_ρ_hist(r::AbstractVector{T}; weights=nothing) where T <: Real
    r_bins = round(Int64, 0.1 * sqrt(length(r)))
    return calc_ρ_hist(r, r_bins, weights=weights)
end


function calc_ρ_hist(snap::Snapshot, bins; weights=snap.masses, x_cen=snap.x_cen)
    r = calc_r(snap, x_cen)
    return calc_ρ_hist(r, bins; weights=weights)
end

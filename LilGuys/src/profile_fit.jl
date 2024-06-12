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
function calc_ρ_hist(r, bins::AbstractVector; weights=nothing)
    if weights == nothing
        weights = ones(length(r))
    end

    counts = Arya.histogram(r, bins, weights=weights, normalization=:count).values

    Vs = 4/3 * π * diff(bins .^ 3)
    return bins, counts ./ Vs
end


function calc_ρ_hist(r, bins::Int; weights=nothing)
    x1 = log10(minimum(r))
    x2 = log10(maximum(r))
    x = LinRange(x1, x2, bins)
    bins = 10 .^ x
    bins = percentile(radii, LinRange(0, 100, Nr+1))
    return calc_ρ_hist(r, bins; weights=weights)
end


function calc_ρ_hist(r; weights=nothing)
    r_bins = round(Int64, 0.1 * sqrt(length(r)))
    return calc_ρ_hist(r, r_bins, weights=weights)
end


function calc_ρ_hist(snap::Snapshot, bins; weights=snap.masses, x_cen=snap.x_cen)
    r = calc_r(snap, x_cen)
    return calc_ρ_hist(r, bins; weights=weights)
end

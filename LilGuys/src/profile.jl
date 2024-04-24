
"""

sorts a snapshot by radius from 0
"""
function sort_by_r(snap::Snapshot)
    return snap[sortperm(calc_r(snap.positions))]
end


"""
Calculates the mass inside each bin
"""
function calc_m(r, masses, r_bins)
    N = length(r_bins) - 1

    ms = zeros(N)
    for i in 1:N
        r_min = r_bins[i]
        r_max = r_bins[i+1]
        ms[i] = sum(masses[r_min .< r .< r_max])
    end

    return ms
end


function calc_ρ_hist(r, bins::AbstractVector; weights=nothing)
    if weights == nothing
        weights = ones(length(r))
    end

    _, counts = calc_histogram(r, bins, weights=weights)

    Vs = 4/3 * π * diff(bins .^ 3)
    return bins, counts ./ Vs
end


"""
    calc_̢ρ_hist(r, r_bins::Int; weights=nothing)

Calculates the density profile given the masses and the radii for a spherical distribution"""
function calc_ρ_hist(r, bins::Int; weights=nothing)
    x1 = log10(minimum(r))
    x2 = log10(maximum(r))
    x = LinRange(x1, x2, bins)
    bins = 10 .^ x
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

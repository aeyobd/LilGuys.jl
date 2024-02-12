using Base: @kwdef
using StatsBase: fit, Histogram, Weights


# check these!
ρ_s(r, r_s, n) = exp(-(r/r_s)^n)
ρ_s_int(r, r_s, n) = -r_s^n * exp(-(r/r_s)^n) * (r/r_s)^n * (n+1)


"""
A 1-d representation of the profile of a galaxy.
"""
@kwdef struct Profile
    snap::Snapshot
    rs::OptVector = nothing
    Ms::OptVector = nothing 
    ρs::OptVector = nothing
    Φs::OptVector = nothing
    Vs_circ::OptVector = nothing
end


function Profile(snap::Snapshot)
    x_c = centroid(snap.positions)
    v_c = centroid(snap.velocities)

    return Profile(snap, x_c, v_c)
end

"""
Given a snaphot, returns a profile object
"""
function Profile(snap::Snapshot, x_c, v_c, Nr=20, Ne=20)
    sorted = copy(snap)
    sorted.positions .-= x_c
    sorted.velocities .-= v_c
    sorted = sort_by_r(sorted)

    N = length(sorted)
    r = calc_r(sorted)
    r_bins = make_equal_number_bins(r, Nr)
    r_mids = (r_bins[1:end-1] + r_bins[2:end]) / 2

    m = calc_m(r, sorted.masses, r_bins)
    M = cumsum(m)

    V_circ = [get_V_circ(sorted, r) for r in r_mids]
    ρs = calc_ρ_profile(r, m, r_bins)

    return Profile(snap=snap, rs=r_mids, Ms=M, 
                   ρs=ρs, Vs_circ=V_circ)
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


function calc_ρ_profile(r, ms, r_bins)
    N = length(r_bins) - 1
    ρs = zeros(N)

    for i in 1:N
        r_min = r_bins[i]
        r_max = r_bins[i+1]
        V = 4/3 * π * (r_max^3 - r_min^3)
        ρs[i] = ms[i] / V
    end

    return ρs
end





"""
calculates equal number bins over the array x with n values per bin.
"""
function make_equal_number_bins(x, n)
    N = length(x)
    xs = sort(x)
    Nbins = round(Int, N/n) + 1
    n = N / Nbins

    bins = zeros(Nbins)

    dx1 = xs[2] - xs[1]
    dx2 = xs[end] - xs[end-1]
    bins[1] = xs[1] - dx1/2
    bins[end] = xs[end] + dx2/2

    for i in 2:(Nbins-1)
        ii = i * n
        xl = xs[floor(Int, ii)]
        xh = xs[ceil(Int, ii)]
        bins[i] = (xl + xh) / 2
    end

    return bins
end



function get_V_circ_max(snap::Snapshot)
    r = calc_r(snap.positions)
    N = length(snap)

    return max(get_V_circ.(snap, r))
end


function get_V_circ(snap::Snapshot, r)
    rs = calc_r(snap)
    M = sum((rs .< r) .* snap.masses)
    return calc_V_circ(M, r)
end



function sort_by_r(snap::Snapshot)
    return snap[sortperm(calc_r(snap.positions))]
end

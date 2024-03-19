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


"""
    Profile(snap)
Given a snapshot `snap`, returns a profile object.
"""
function Profile(snap::Snapshot)
    x_c = centroid(snap.positions)
    v_c = centroid(snap.velocities)

    return Profile(snap, x_c, v_c)
end


function write(file::IO, profile::Profile)
    write(file, profile.rs)
    write(file, profile.Ms)
    write(file, profile.ρs)
    write(file, profile.Φs)
    write(file, profile.Vs_circ)
end

"""
    Profile(snap, x_c, v_c; Nr=20, Ne=20)

Given a snaphot, returns a profile object
# Arguments
- `snap::Snapshot`: the snapshot to be profiled
"""
function Profile(snap::Snapshot, x_c, v_c; Nr=20, Ne=20)
    sorted = copy(snap)
    sorted.positions .-= x_c
    sorted.velocities .-= v_c
    sorted = sort_by_r(sorted)
    filt = calc_E_spec(sorted) .< 0
    sorted = sorted[filt]

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


"""Returns a centred snapshot"""
function centre_obsolete(snap::Snapshot)
    x_c = centroid(snap.positions)
    v_c = centroid(snap.velocities)

    centred = copy(snap)
    centred.positions .-= x_c
    centred.velocities .-= v_c

    return centred
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


"""Calculates the densities given the masses and the radii for a spherical distribution"""
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




"""Calculates the maximum circular velocity of a snapshot"""
function get_V_circ_max(snap::Snapshot)
    r = calc_r(snap.positions)
    N = length(snap)

    return max(get_V_circ.(snap, r))
end


"""calculates the circular velocity for each radius of a particle in snap"""
function get_V_circ(snap::Snapshot, r)
    rs = calc_r(snap)
    M = sum((rs .< r) .* snap.masses)
    return calc_V_circ(M, r)
end


"""sorts a snapshot by radius from 0"""
function sort_by_r(snap::Snapshot)
    return snap[sortperm(calc_r(snap.positions))]
end

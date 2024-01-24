using Base: @kwdef

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


function Profile(snap::Snapshot, r_s, n, Nr=20, Ne=20)
    sorted = sort_by_r(center(snap))

    N = length(sorted)
    r = calc_r(sorted)
    r_bins = create_r_bins(r[1], r[-1], Nr)
    r_mids = (r_bins[1:end-1] + r_bins[2:end]) / 2

    m = snap.masses
    M = collect(1:N) .* m
    ν = get_ν(r, m, r_bins)
    V_circ = map(r->get_V_circ(sorted, r), r_bins)
    ρs = calc_ρ_profile(r_s, M, r_bins)

    return Profile(snap=snapshot, rs=r_mids, Ms=M, 
                   ρs=ρs, V_circ=V_circ)
end


function calc_ρ_profile(r_s, ms, r_bins)
    N = length(r_bins) - 1

    ρs = zeros(N)
    for i in 1:N
        r_min = r_bins[i]
        r_max = r_bins[i+1]
        M = sum(ms[r_min .< rs .< r_max])
        V = 4/3 * π * (r_max^3 - r_min^3)
        ρs[i] = M / V
    end

    return ρs
end

function centre(snap::Snapshot)
    cen = ss_center(snap)
    snap1 = copy(snap)
    snap1.postions .-= cen.x_c
    snap1.velocities .-= cen.v_c
    return snap1
end


function create_r_bins(r_min, r_max, n_bins)
    return exp10.(range(log10(r_min), stop=log10(r_max), length=n_bins+1))
end


function get_ν(r, m, r_bins)
    hist = fit(Histogram, r, m, r_bins)
    V = r_bins[2:end] .^ 3 - r_bins[1:end-1] .^ 3
    V .*= 4/3 * π
    return hist.weights * m ./ V 
end


# check these!
ρ_s(r, r_s, n) = exp(-(r/r_s)^n)
ρ_s_int(r, r_s, n) = -r_s^n * exp(-(r/r_s)^n) * (r/r_s)^n * (n+1)



function get_V_circ_max(snap::Snapshot)
    r = calc_r(snap.pos)
    M = arange(length(snapshot)) .* snapshot.m
    return max(get_V_circ(M, r))
end

function get_V_circ(snap::Snapshot, r)
    rs = calc_r(snap)
    M = sum(rs < r) * snap.m
    return V_circ(M, r)
end



function sort_by_r(snap::Snapshot)
    return snap[sortperm(calc_r(snap.pos))]
end

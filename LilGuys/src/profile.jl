using Base: @kwdef

"""
A 1-d representation of the profile of a galaxy.
"""
@kwdef struct Profile
    cen::FuzzyPhase
    r::OptVector = nothing
    M::OptVector = nothing
    ρ::OptVector = nothing
    V_circ::OptVector = nothing
end


function Profile(snap::Snapshot, r_s, n, Nr=20, Ne=20)
    sorted = sort_by_r(center(snap))

    N = length(sorted)
    r = get_r(sorted)
    r_bins = create_r_bins(r[1], r[-1], Nr)

    m = snap.m
    M = collect(1:N) .* m
    ν = get_ν(r, m, r_bins)
    V_circ = map(r->get_V_circ(sorted, r), r_bins)
    return Profile(r, M, ν, V_circ)
end


function Base.getproperty(p::Profile, s::Symbol)
    if s ∈ fieldnames(Profile)
        return getfield(p, s)
    elseif s ∈ fieldnames(Snapshot)
        return getfield(p.snap, s)
    else
        error("Profile has no field $s")
    end
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


get_r(snapshot::Snapshot) = get_r(snapshot.pos)

# check these!
ρ_s(r, r_s, n) = exp(-(r/r_s)^n)
ρ_s_int(r, r_s, n) = -r_s^n * exp(-(r/r_s)^n) * (r/r_s)^n * (n+1)




function get_V_circ_max(snap::Snapshot)
    r = get_r(snap.pos)
    M = arange(length(snapshot)) .* snapshot.m
    return max(get_V_circ(M, r))
end

function get_V_circ(snap::Snapshot, r)
    rs = get_r(snap)
    M = sum(rs < r) * snap.m
    return V_circ(M, r)
end



function sort_by_r(snap::Snapshot)
    return snap[sortperm(get_r(snap.pos))]
end

# Finding the center of a snapshot
# This process turns out to be more nuanced, so 
# we have several methods here
#
# Mean
#   The mean center is the most straightforward, but
#   during late-stage tidal stripping, this stops being as useful.
#   The uncertanty is simply the standard error on the mean.
#
# Shrinking Spheres
#
#
# Potential center
#   
#
# Statistical center
#   By combining the above methods, we can find a bayesian estimate
#   of the center given a prior. 
#
#


import StatsBase as sb
import SpecialFunctions: erf
import NearestNeighbors as nn
import Interpolations: linear_interpolation, Line



"""Finds the centre using the shrinking sphere method"""
function ss_centre(snap::Snapshot, perc=95, min_frac=0.2)
    N = length(snap)
    N_min = ceil(min_frac * N)

    snap1 = copy(snap)
    idx0 = argmin(snap1.Φ)
    xs = []
    vs = []


    x0 = snap1.pos[:, idx0]
    v0 = snap1.vel[:, idx0]

    while length(snap1) > N_min
        push!(xs, x0)
        push!(vs, v0)

        snap1 = cut_outside(snap1, x0, perc)
        calc_Φ!(snap1)
        snap1 = cut_unbound(snap1, x0, v0)

        x0, v0 = centroid(snap1.pos)
    end

    print(std(xs))
    print(std(vs))

    return x0, v0
end





function phase_volume(snap::Snapshot; k=5)
    δrs = []
    δvs = []
    tree = nn.KDTree(snap.pos)
    for i in 1:length(snap)
        δr, δv = phase_volume(snap, tree, i)
        push!(δrs, δr)
        push!(δvs, δv)
    end
    return δrs, δvs
end


function phase_volume(snap, tree, idx; k=10, s=5, η=0.3)
    idxs, dists = nn.knn(tree, snap.pos[:, idx], k)
    filt = dists .> 0
    idxs = idxs[filt]
    dists = dists[filt]

    r_std = maximum(dists) / sqrt(k)
    vs = r(snap.vel[:, idxs] .- snap.vel[:, idx])
    v = sb.mean(sort(vs)[1:s])
    v_std = η * v / sqrt(s)
    return r_std, v_std
end




function fuzzy_centre(snap::Snapshot, p::PhasePoint; percen=5, )
    cen = potential_centre(snap, percen=percen)
    δr, δv = phase_volume(snap, p, k=10)

    weights = ones(length(snap))

    for i in 1:10
        weights = bound_probabilities(snap, k=10)
    end
end


"""
given a centered snapshot, returns a interpolated potential a a function 
of r"""
function calc_radial_Φ(snap::Snapshot)
    # work outside in 
    rs = sort(r(snap.pos))
    m = snap.m

    N = length(snap)
    M_inside = m * collect(1:N)

    Φ_inside = -G * M_inside ./ rs
    
    Φ_outside = -G*snap.m .* [ sum(1 ./ rs[i:end]) for i in 1:N+1]

    pushfirst!(rs, 0)
    pushfirst!(Φ_inside, 0)

    Φs = Φ_inside .+ Φ_outside

    return linear_interpolation(rs, Φs, extrapolation_bc=Line())
end





function bound_probabilities(snap::Snapshot; k=5)
    Φ = calc_radial_Φ(snap)
    δxs, δvs = knn_err(snap, k=k)

    probs = zeros(length(snap))

    for i in 1:length(snap)
        pos1 = snap.pos[:, i]
        vel1 = snap.vel[:, i]

        δr = δxs[i]
        δv = δvs[i]
        
        r0 = r(pos1)
        v0 = r(vel1)
        e0 = E_spec(Φ(r0), v0)
        eh = E_spec(Φ(r0 + δr), v0 + δv)
        el = E_spec(Φ(max(r0 - δr, 0)), max(v0 - δv, 0))

        δe = (eh - el) / 2

        # how much of the probability distribution is less than zero
        probs[i] = normal_cdf(0, e0, δe)
    end

    return probs
end



### Utilities for finding the center
#
function centroid(x::Matrix{F}, 
        weights::Vector{F}=ones(length(x))
    return sum(x .* w, dims=2) /sum(weights)
end

function centroid_err(x::Matrix{F}, 
        weights::Vector{F}=ones(length(x))
    variance= sum((x.^2 .* w), dims=2) / sum(weights)
    return sqrt(variance)
end


function centroid(snap::Snapshot, weights::Vector{F}=ones(length(snap))
    pos_c = centroid(snap.pos, weights)
    vel_c = centroid(snap.vel, weights)
    δx = centroid_err(snap.pos .- pos_c, weights)
    δv = centroid_err(snap.vel .- vel_c, weights)

    return FuzzyPhase(pos_c, vel_c, δx, δv)
end


function potential_centre(snap::Snapshot; percen=5)
    threshhold = percentile(snap.Φ, percen)
    filt = snap.Φ .< threshhold
    return centroid(snap[filt])
end



function cut_outside(snap, cerc=95)
    r = r(snap.pos)
    r_cut = percentile(r, perc)
    println(r_cut)
    filt = r .<= r_cut
    return snap[filt]
end


function cut_unbound(snap::Snapshot, x0, v0)
    x1 = snap.pos .- x0
    v1 = snap.vel .- v0

    E_kin = 0.5 * r(v1).^2
    E_spec = snap.Φ .+ E_kin

    filt = E_spec .< 0
    return snap[filt]
end


function normal_cdf(x, μ, σ)
    z = (x - μ) / σ
    return 0.5 * (1 + erf(z / sqrt(2)))
end



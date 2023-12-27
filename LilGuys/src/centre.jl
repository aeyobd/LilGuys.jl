import StatsBase as sb
import SpecialFunctions: erf
import NearestNeighbors as nn


"""
Computes the centre of a snapshot while accounting for uncertanties

Parameters
----------
snap : Snapshot
    The snapshot to find the centre of
threshold : Float
    The threshold for the probability distribution. Any point with a 
    probability less than this will be ignored.
γ : Float
    Momentum decay on centres.
"""
function fuzzy_centre!(snap_i::FSnapshot; min_fraction=0.1, threshold=0.1, 
        β=0.9, dx_min=0.05, dv_min=0.001, max_iter=10)
    snap = copy(snap_i)
    dcen = copy(cen)

    for _ in 1:max_iter
        update_weights!(snap, threshold=threshold, β=β)
        cen = centroid(snap_c, snap_c.m .* snap_c.w)

        dx = β*dx + (1-β) * cen.pos
        dv = β*dv + (1-β) * cen.vel
        δx = β*dcen.δx + (1-β) * (norm(dx) + cen.δx)
        δv = β*dcen.δv + (1-β) * (norm(dv) + cen.δv)
        dcen = FuzzyPhase(dx, dv, δx, δv)

        update_centre!(snap_c, dcen)
        frac = mean(snap_c.w .> threshold)
        xc = xc .+ dx
        vc = vc .+ dv

        if (norm(dx) < dx_min) && (norm(dv) < dv_min)
            break
        end

        println("f = ", frac)
        println(cen)
        if frac < min_fraction
            break
        end
    end

    cen = FuzzyPhase(xc, vc, cen.δx, cen.δv)
    return cen, snap_c.w
end


function update_weights!(snap::FSnapshot; β=0.5, threshold=0)
    weights = bound_probabilities(snap, k=10)
    weights = β * snap.w .+ (1-β) * weights
    weights[weights .< threshold] .= 0
    snap.w = weights
    return snap
end


function update_centre!(snap::FSnapshot, p::FuzzyPhase; percen=0.5)
    cen = potential_centre(snap, percen=percen)
    δrs, δvs = phase_volumes(snap, k=10)
    δrs .= cen.δx
    δvs .+= cen.δv
    snap.pos .-= cen.pos 
    snap.vel .-= cen.vel
    snap.δr = δrs
    snap.δv = δvs
    return snap
end




function bound_probabilities(snap::Snapshot; k=5)
    Φ = calc_radial_Φ(snap)
    δxs, δvs = phase_volumes(snap, k=k)

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
function centroid(x::Matrix{F}, weights::Vector{F})
    w = reshape(weights, :, 1) ./ sum(weights)
    return (x * w)[:, 1]
end

function centroid(x::Matrix{F})
    c =  sum(x, dims=2) / size(x, 2)
    return c[:, 1]
end

function centroid_err(x::Matrix{F}, weights::Vector{F})
    w = reshape(weights, :, 1) ./ sum(weights)
    variance = mean(x.^2 * w)
    return sqrt(variance)
end

function centroid_err(x::Matrix{F})
    variance= mean(sum(x.^2, dims=2) / size(x, 2))
    return sqrt(variance)
end


function centroid(snap::AbstractSnapshot)
    pos_c = centroid(snap.pos)
    vel_c = centroid(snap.vel)
    δr = centroid_err(snap.pos .- pos_c)
    δv = centroid_err(snap.vel .- vel_c)

    return FuzzyPhase(pos_c, vel_c, δr, δv)
end

function centroid(snap::AbstractSnapshot, weights::Vector{F})
    pos_c = centroid(snap.pos, weights)
    vel_c = centroid(snap.vel, weights)
    δr = centroid_err(snap.pos .- pos_c, weights)
    δv = centroid_err(snap.vel .- vel_c, weights)

    return FuzzyPhase(pos_c, vel_c, δr, δv)
end


function potential_centre(snap::AbstractSnapshot; percen=5)
    threshhold = percentile(snap.Φ, percen)
    filt = snap.Φ .< threshhold
    return centroid(snap[filt])
end



function normal_cdf(x, μ, σ)
    z = (x - μ) / σ
    return normal_cdf(z)
end

function normal_cdf(z)
    return 0.5 * (1 + erf(z / sqrt(2)))
end



function phase_volumes(snap::AbstractSnapshot; k=5)
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


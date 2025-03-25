# A work in progress method for bayesian centre finding for the 
# time series of snapshots, so especially for outputs.

import NearestNeighbors as nn
import SpecialFunctions: gamma

import Base: @kwdef


@kwdef struct FuzzyCentreState
    weights::Vector{F} # not sure...  binom?
    δr::OptVector = nothing #
    δv::OptVector = nothing
    w::OptVector = nothing
end




"""
Computes the centre of a snapshot while accounting for uncertanties.

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
function fuzzy_centre!(snap_i::Snapshot; min_fraction=0.1, threshold=0.1, 
        β=0.9, dx_min=0.05, dv_min=0.001, max_iter=10)
    snap = copy(snap_i)
    dcen = copy(cen)

    for _ in 1:max_iter
        update_weights!(snap, threshold=threshold, β=β)
        cen = centroid(snap_c, snap_c.m .* snap_c.w)

        dx = β*dx + (1-β) * cen.positions
        dv = β*dv + (1-β) * cen.velocity
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


function update_weights!(snap::Snapshot; β=0.5, threshold=0)
    weights = bound_probabilities(snap, k=10)
    weights = β * snap.w .+ (1-β) * weights
    weights[weights .< threshold] .= 0
    snap.w = weights
    return snap
end


function update_centre!(snap::Snapshot, p; percen=0.5)
    cen = potential_centre(snap, percen=percen)
    δrs, δvs = phase_volumes(snap, k=10)
    δrs .= cen.δx
    δvs .+= cen.δv
    snap.positions .-= cen.pos 
    snap.velocities .-= cen.vel
    snap.δr = δrs
    snap.δv = δvs
    return snap
end




function bound_probabilities(snap::Snapshot; k=5)
    Φ = potential_radial_func(snap)
    δxs, δvs = phase_volumes(snap, k=k)

    probs = zeros(length(snap))

    for i in 1:length(snap)
        pos1 = snap.positions[:, i]
        vel1 = snap.velocities[:, i]

        δr = δxs[i]
        δv = δvs[i]
        
        r0 = radii(pos1)
        v0 = radii(vel1)
        e0 = specific_energy(Φ(r0), v0)
        eh = specific_energy(Φ(r0 + δr), v0 + δv)
        el = specific_energy(Φ(max(r0 - δr, 0)), max(v0 - δv, 0))

        δe = (eh - el) / 2

        # how much of the probability distribution is less than zero
        probs[i] = normal_cdf(0, e0, δe)
    end

    return probs
end


"""
Computes the sizes of each particle in phase space.
I.e. the 

Parameters
----------
k : Int
    The number of nearest neighbours to use
s : Int
    The number of nearest neighbours to use for the velocity
η : Float
    The fraction of the velocity to use for the velocity standard deviation
"""
function phase_volumes(snap::Snapshot; k=5, kwargs...)
    tree = nn.KDTree(snap.positions)
    idxs, dists = nn.knn(tree, snap.positions, k+1, true)

    δrs = []
    δvs = []

    gamma_ratio = gamma(k+1/3) / gamma(k) 

    for i in 1:length(snap)
        idx = idxs[i][2:end]
        rs = dists[i][2:end]

        vel_0 = snap.velocities[:, i]
        vs = radii(snap.velocities[:, idx], vel_0)

        δr, δv = _phase_volume(rs, vs; gamma_ratio=gamma_ratio, kwargs...)
        push!(δrs, δr)
        push!(δvs, δv)
    end
    return δrs, δvs
end


@doc raw"""
    _phase_volume(rs::Vector{F}, vs::Vector{F}; η=1.0, s=nothing, gamma_ratio=nothing)

Returns the 1D velocity dispersion and the characteristic nearest-neighbor radius of the particle given positions and velocities for the kth nearest neighbors (sorted by position). 

Given a poisson process, the expected distance to the nth nearest neighbour in 3D is
``math
r_n = \frac{\Gamma(k+1/3)}{\Gamma(k)} \left( \frac{3}{4\pi\rho} \right)^{1/3}
``
where \(k\) is the number of nearest neighbours, and \(\rho\) is the density of the particles.
As such, we return the value
``math
r_m = r_k \left( \frac{4\pi}{3} \right)^{1/3} \frac{\Gamma(k)}{\Gamma(k+1/3)}
``
such that $\frac{1}{4\pi/3 r_m^3}$ is a (probably biased) estimate of the local density. 
"""
function _phase_volume(rs::Vector{F}, vs::Vector{F}; η=1.0, s=nothing, gamma_ratio=nothing)
    k = length(rs)

    if s === nothing
        s = k
    end

    r_std = rs[end] / gamma_ratio

    v2 = sort(vs)[1:s] .^ 2
    v_std = η * sqrt(sum(v2) / s)

    return r_std, v_std
end


import NearestNeighbors as nn
import Optim

function centre_ρ(snap::Snapshot; k=5)
    tree = nn.KDTree(snap.positions)
    ρs = calc_ρs(snap; k=k)
    p0 = snap.positions[:, argmax(ρs)]

    f(p) = -calc_ρ(p, tree, snap.masses; k=k) # find maximum
    optimizer = Optim.LBFGS()
    result = Optim.optimize(f, p0, optimizer)
    return result
end


"""
Computes the densities at each particle in a snapshot

"""
function calc_ρs(positions::AbstractMatrix, masses::AbstractVector; k=5, kwargs...)

    tree = nn.KDTree(positions)
    idxs, dists = nn.knn(tree, positions, k, true)
    N = size(positions, 2)

    ρs = Vector{F}(undef, N)

    for i in 1:N
        rs = dists[i]
        idx = idxs[i]
        ms = masses[idx]

        ρs[i] = _calc_ρ(rs, ms; kwargs...)
    end

    return ρs
end


function calc_ρs(snap::Snapshot; kwargs...)
    return calc_ρs(snap.positions, snap.masses; kwargs...)
end


function calc_ρ(position, tree::nn.KDTree, masses::AbstractVector; k=5, kwargs...)

    idxs, dists = nn.knn(tree, position, k, true)

    ms = masses[idxs]

    ρ = _calc_ρ(dists, ms; kwargs...)

    return ρ
end



function _calc_ρ(rs::Vector{F}, ms=ones(length(rs)); η=0.3, h=nothing, s=nothing)
    if h === nothing
        h = η * rs[end]
    end

    ws = sph_kernel.(rs, h)
    ρ = sum(ws .* ms)
    return ρ
end


function sph_kernel(q)
    σ =  21 / (16π)

    if 0 <= q < 2
        return σ *  (1 - q/2)^4 * (2q + 1)
    else
        return 0
    end
end


function sph_kernel(r, h)
    q = r / h
    return 1/h^3 * sph_kernel(q)
end

Point = Vector{F}
abstract type AbstractState end

struct Centre
    position::Point
    position_err::F
    velocity::Point
    velocity_err::F
    acceleration::Point
    acceleration_err::F
end

@kwdef mutable struct StaticState <: AbstractState
    centre::Centre
    method = "com"
end


function StaticState(snap::Snapshot; method="com")
    return StaticState(centre=Centre(), method=method)
end

function calc_centre!(state::StaticState, snap)
    if state.method == "com"
        state.centre = centre_of_mass(snap)
    elseif state.method == "potential"
        state.centre = centre_weighted_potential(snap)
    end
    return state
end


function calc_next_centre!(state::StaticState, snap)
    calc_centre!(state, snap)
end


function Centre()
    return Centre(zeros(3), NaN,
                  zeros(3), NaN,
                  zeros(3), NaN)
end




function weighted_centre(snap::Snapshot, weights::AbstractVector)
    position = centroid(snap.positions, weights)
    position_err = centroid_err(snap.positions, weights)
    velocity = centroid(snap.velocities, weights)
    velocity_err = centroid_err(snap.velocities, weights)
    acceleration = centroid(snap.accelerations, weights)
    acceleration_err = centroid_err(snap.accelerations, weights)

    return Centre(position, position_err, velocity, velocity_err, acceleration, acceleration_err)
end


function centre_of_mass(snap::Snapshot)
    return weighted_centre(snap, snap.masses)
end


function centre_weighted_potential(snap::Snapshot)
    return weighted_centre(snap, -snap.Φs)
end


function centre_potential_percen(snap::Snapshot, percen=5)
    Φcut = percentile(snap.Φs, percen)
    filt = snap.Φs .< Φcut
    return centre_of_mass(snap[filt])
end
    


function ϕ_eff(positions, velocities, masses; h=0.08)
	N = length(masses)
	ϕ = zeros(N)

	for i in 1:N
		if i % 100 == 0
			print("\r $i / $N")
		end

		ϕ_g = @. -masses / sqrt(calc_r(positions[:, i], positions)^2 + h^2)
		ϕ_v =  1/N * 1/2 * calc_r(velocities[:, i], velocities).^2
			
		ϕ[i] = sum(min.(ϕ_g .+ ϕ_v, 0))
		
	end
	return ϕ

end

function calc_ρ_eff(x_vec, v_vec, positions, velocities, masses; h)
    N = length(masses)

    ϕ_g = @. -masses / sqrt(calc_r(x_vec, positions)^2 + h^2)
    ϕ_v =  1/N * 1/2 * calc_r(v_vec, velocities).^2
    ϕs = @. min(ϕ_g + ϕ_v, 0)
    return sum(ϕs)
end





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


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


function mean_centre(snap::Snapshot, filter)
    position = centroid(snap.positions[:, filter])
    position_err = centroid_err(snap.positions[:, filter])
    velocity = centroid(snap.velocities[:, filter])
    velocity_err = centroid_err(snap.velocities[:, filter])
    if snap.accelerations !== nothing
        acceleration = centroid(snap.accelerations[:, filter])
        acceleration_err = centroid_err(snap.accelerations[:, filter])
    else
        acceleration = zeros(3)
        acceleration_err = NaN
    end

    return Centre(position, position_err, velocity, velocity_err, acceleration, acceleration_err)
end





function weighted_centre(snap::Snapshot, weights::AbstractVector)
    position = centroid(snap.positions, weights)
    position_err = centroid_err(snap.positions, weights)
    velocity = centroid(snap.velocities, weights)
    velocity_err = centroid_err(snap.velocities, weights)
    if snap.accelerations !== nothing
        acceleration = centroid(snap.accelerations, weights)
        acceleration_err = centroid_err(snap.accelerations, weights)
    else
        acceleration = zeros(3)
        acceleration_err = NaN
    end

    return Centre(position, position_err, velocity, velocity_err, acceleration, acceleration_err)
end


function centre_of_mass(snap::Snapshot)
    return weighted_centre(snap, snap.masses)
end


function centre_weighted_potential(snap::Snapshot)
    return weighted_centre(snap, -snap.Φs)
end


function centre_potential_percen(snap::Snapshot, percen=5)
    if snap.Φs === nothing
        Φs = calc_radial_discrete_Φ(snap)
    else
        Φs = snap.Φs
    end
    Φcut = percentile(Φs, percen)
    filt = Φs .< Φcut
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






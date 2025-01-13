
Point = Vector{F}

"""Abstract centre finding state"""
abstract type AbstractState end

"""
A struct representing a centre estimate for a snapshot
"""
struct Centre
    position::Point
    position_err::F
    velocity::Point
    velocity_err::F
    acceleration::Point
    acceleration_err::F
end


function Centre()
    return Centre(zeros(3), NaN,
                  zeros(3), NaN,
                  zeros(3), NaN)
end



"""State representing a static centre"""
@kwdef mutable struct StaticState <: AbstractState
    centre::Centre
    method = "com"
end


function StaticState(snap::Snapshot; method="com", verbose=false)
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







function centre_of_mass(snap::Snapshot)
    return weighted_centre(snap, snap.masses)
end


function centre_weighted_potential(snap::Snapshot)
    return weighted_centre(snap, -snap.Φs)
end


"""
Simple weighted centre of a snapshot
"""
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



"""
    centre_potential_percen(snap, q)

Computes the COM of particles with a potential below a certain quantile q.
"""
function centre_potential_percen(snap::Snapshot, q=0.05)
    if snap.Φs === nothing
        Φs = calc_radial_discrete_Φ(snap)
    else
        Φs = snap.Φs
    end
    Φcut = quantile(Φs, q)
    filt = Φs .< Φcut
    return centre_of_mass(snap[filt])
end


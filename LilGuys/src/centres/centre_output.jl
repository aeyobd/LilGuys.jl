"""
    calc_centres(StateType, out::Output; kwargs...)

Calculates the centres for each snapshot in the output.
The details of the implementation are based on the given StateType.
"""
function calc_centres(StateType, out::Output; reinit_state=false, skip=1, kwargs...)
    state = StateType(out[1]; kwargs...)

    calc_centre!(state, out[1])
    centres = [state.centre]

    time = out.times[1]

    for i in 2:skip:length(out)
        dt = out.times[i] - time

        update_prior!(state, dt)
        cen = state.centre
        if reinit_state
            state = StateType(out[i]; kwargs...)
            state.centre = cen
            calc_centre!(state, out[i])
        else
            calc_next_centre!(state, out[i])
        end

       
        push!(centres, state.centre)
        time = out.times[1]
    end

    return centres
end



"""
    calc_centre(StateType, snap::Snapshot; kwargs...)

Calculates the centres of a snapshot.
The details of the implementation are based on the given StateType.
"""
function calc_centre(StateType, snap::Snapshot; kwargs...)
    state = StateType(snap; kwargs...)
    calc_centre!(state, snap)
    return state.centre
end



function update_prior!(state::AbstractState, dt::Real)
    state.centre = leapfrog(state.centre, dt)
end


function leapfrog(centre, dt)
    position = centre.position .+ centre.velocity * dt / 2
    velocity = centre.velocity .+ centre.acceleration * dt
    position .+= centre.velocity * dt / 2

    position_err = centre.position_err + centre.velocity_err * dt / 2
    velocity_err = centre.velocity_err + centre.acceleration_err * dt
    position_err += centre.velocity_err * dt / 2

    return Centre(position, position_err, velocity, velocity_err, centre.acceleration, centre.acceleration_err)
end



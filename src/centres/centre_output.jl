"""
    calc_centres(StateType, out::Output; kwargs...)

Calculates the centres for each snapshot in the output.
The details of the implementation are based on the given StateType.
"""
function calc_centres(StateType, out::Output; reinit_state=false, verbose=false, skip=1, kwargs...)

    out_idx = 1:skip:length(out)

    centres = Vector{Centre}(undef, length(out_idx))

    state = StateType(out[1]; verbose=verbose, kwargs...)

    calc_centre!(state, out[1])
    centres[1] = state.centre

    time_last = out.times[1]

    for ii in 2:length(out_idx)
        i = out_idx[ii]
        i_last = out_idx[ii-1]
        dt = out.times[i] - out.times[i_last]
        snap = out[i]

        update_prior!(state, dt)
        cen = state.centre

        if verbose
            println("using prior: $(cen.position)")
            println("error : $(cen.position_err)")
        end

        if reinit_state
            state = StateType(snap; verbose=verbose, kwargs...)
            state.centre = cen
            calc_centre!(state, snap)
        else
            calc_next_centre!(state, snap)
        end

       
        centres[ii] = state.centre

        if verbose
            println("completed snapshot $i")
        end
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



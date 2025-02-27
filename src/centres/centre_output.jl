"""
    calc_centres(StateType, out::Output; reinit_state=false, skip=1, kwargs...)

Calculates the centres for each snapshot in the output.
The details of the implementation are based on the given StateType.
Reinit_state causes a new instance of the state to be reinitialised for each snapshot (instead of using information from previous snapshot). Skip sets the interval between snapshots to calculate the centres.
"""
function calc_centres(StateType::Type, out::Output; reinit_state=false, skip=1, kwargs...)

    out_idx = 1:skip:length(out)

    centres = Vector{Centre}(undef, length(out_idx))

    state = StateType(out[1]; kwargs...)

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

        @info "using prior: $(cen.position)"
        @info "error : $(cen.position_err)"

        if reinit_state
            state = StateType(snap; kwargs...)
            state.centre = cen
            calc_centre!(state, snap)
        else
            calc_next_centre!(state, snap)
        end

       
        centres[ii] = state.centre

        @info "completed snapshot $i"
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


"""
    update_prior!(state, dt::Real)

Updates the state's centre using the leapfrog method.
"""
function update_prior!(state::AbstractState, dt::Real)
    state.centre = leapfrog(state.centre, dt)
end


"""
    leapfrog(centre, dt)

Updates the centres position and velocity using the leapfrog method.
"""
function leapfrog(centre::Centre, dt::Real)::Centre
    position = centre.position .+ centre.velocity * dt / 2
    velocity = centre.velocity .+ centre.acceleration * dt
    position .+= centre.velocity * dt / 2

    position_err = centre.position_err + centre.velocity_err * dt / 2
    velocity_err = centre.velocity_err + centre.acceleration_err * dt
    position_err += centre.velocity_err * dt / 2

    return Centre(position, position_err, velocity, velocity_err, centre.acceleration, centre.acceleration_err)
end



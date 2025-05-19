import Base: -, reverse
import CSV

"""
Represents an Orbit. All quantities are in code units.
"""
Base.@kwdef struct Orbit
    "Time in code units"
    times::Vector{F}

    "Position in code units"
    positions::Matrix{F}

    "Velocity in code units"
    velocities::Matrix{F}

    "Acceleration"
    accelerations::Union{Matrix{F}, Nothing} = nothing

    "pericentre with respect to origin"
    pericenter::F = minimum(radii(positions))

    "apocentre with respect to origin"
    apocenter::F = maximum(radii(positions))

    function Orbit(times, positions, velocities, accelerations, pericenter, apocenter)
        @assert_3vector positions
        @assert_3vector velocities
        @assert_same_size positions velocities 

        if length(times) != size(positions, 2) 
            throw(ArgumentError("times and positions should have same number of entries"))
        end

        if !isnothing(accelerations)
            @assert_3vector accelerations
            @assert_same_size positions accelerations
        end


        return new(times, positions, velocities, accelerations, pericenter, apocenter)
    end
end

positions(a::Orbit) = a.positions
velocities(a::Orbit) = a.velocities
accelerations(a::Orbit) = a.accelerations
times(a::Orbit) = a.times

pericenter(a::Orbit) = a.pericenter
apocenter(a::Orbit) = a.apocenter
Base.length(a::Orbit) = length(a.times)


""" Loads an orbit from a CSV file. Expected columns are `time`, `x`, `y`, `z`, `v_x`, `v_y`, `v_z`. """
function Orbit(filename::String)
    df = CSV.read(filename, DataFrame)
    positions = hcat(df.x, df.y, df.z)'
    velocities = hcat(df.v_x, df.v_y, df.v_z)'
    if "a_x" ∈ names(df)
        accelerations =  hcat(df.a_x, df.a_y, df.a_z)'
    else
        accelerations = nothing
    end

    return Orbit(
        times = df.t,
        positions = positions,
        velocities = velocities,
        accelerations = accelerations,
    )
end


""" Converts an orbit to a DataFrame """
function to_frame(a::Orbit)
    df = DataFrame(
        t = a.times,
        x = a.positions[1, :],
        y = a.positions[2, :],
        z = a.positions[3, :],
        v_x = a.velocities[1, :],
        v_y = a.velocities[2, :],
        v_z = a.velocities[3, :]
    )

    if a.accelerations !== nothing
        df = hcat(df, DataFrame(
            a_x = a.accelerations[1, :],
            a_y = a.accelerations[2, :],
            a_z = a.accelerations[3, :]
        ))
    end

    return df
end


""" Writes an orbit to a CSV file """
function Base.write(filename::String, a::Orbit)
    df = to_frame(a)
    CSV.write(filename, df)
end


""" Returns a new orbit which is `a` from `b`'s perspective (`b` is at 0,0). Note acceleration/etc. not currently implemented."""
function (-)(a::Orbit, b::Orbit)
    @assert isapprox(a.times, b.times, rtol=1e-6)

    return Orbit(times=a.times, positions=a.positions - b.positions, velocities=a.velocities - b.velocities)
end


"""
    reverse(a::Orbit)

Reverses the order in time of an orbit
"""
function reverse(a::Orbit)
    time = reverse(a.times)
    positions = reverse(a.positions, dims=2)
    velocities = reverse(a.velocities, dims=2)
    accelerations = a.accelerations === nothing ? nothing : reverse(a.accelerations, dims=2)
    return Orbit(times=time, positions=positions, velocities=velocities, accelerations=accelerations)
end


"""
    resample(a::Orbit, time::AbstractVector)
Resamples an orbit to a new time array
"""
function resample(a::Orbit, time::AbstractVector{<:Real})
    if a.time[2] < a.time[1]
        a = reverse(a)
    end
    x = LilGuys.lerp(a.time, a.positions[1, :]).(time)
    y = LilGuys.lerp(a.time, a.positions[2, :]).(time)
    z = LilGuys.lerp(a.time, a.positions[3, :]).(time)
    v_x = LilGuys.lerp(a.time, a.velocities[1, :]).(time)
    v_y = LilGuys.lerp(a.time, a.velocities[2, :]).(time)
    v_z = LilGuys.lerp(a.time, a.velocities[3, :]).(time)

    return Orbit(time=time, positions=[x y z]', velocities=[v_x v_y v_z]')
end


"""
    leap_frog(gc, acceleration; params...)

Computes the orbit of a Galactocentric object using the leap frog method.

Parameters:
- `gc::Galactocentric`: the initial conditions
- `acceleration::Function`: the acceleration function. should take a position vector, velocity vector, and the current time (for maximum flexibility).
- `dt_max::Real=0.1`: the maximum timestep
- `dt_min::Real=0.001`: the minimum timestep, will exit if timestep is below this
- `time::Real=-10/T2GYR`: the time to integrate to
- `timebegin::Real=0`: the time to start integrating from
- `timestep::Symbol=:adaptive`: the timestep to use, either `:adaptive` or a real number
- `η::Real=0.01`: the adaptive timestep parameter. Adaptive timestep is based on sqrt(η/a) where a is the magnitude of the current acceleration.

"""
function leap_frog(gc, f_acc; 
        dt_max=0.1, dt_min=0.001, timebegin=0, time=-10/T2GYR, 
        timestep=:adaptive, η=0.01, reuse_acceleration=true
    )

    if timestep isa Real
        dt_min = timestep
    end

    t = timebegin

    Nt = round(Int, abs(time / dt_min))

    # setup matricies
    positions = Vector{Vector{Float64}}()
    velocities = Vector{Vector{Float64}}()
    accelerations = Vector{Vector{Float64}}()
    times = Float64[]

    push!(positions, [gc.x, gc.y, gc.z])
    push!(velocities, [gc.v_x, gc.v_y, gc.v_z] / V2KMS)
    push!(accelerations, f_acc(positions[1], velocities[1], t))
    push!(times, 0.)
    is_done = false
    backwards = time < 0

    for i in 1:Nt
        pos = positions[i]
        vel = velocities[i]
        acc = accelerations[i]
        
        if timestep == :adaptive
            dt = min(sqrt(η / radii(acc)), dt_max)
        elseif timestep isa Real
            dt = timestep
        end
        
        if backwards
            dt *= -1
        end
        if abs(dt) < dt_min
            @warn "timestep below minimum timestep"
            is_done = true
        end

        # last step ends at end of range
        if (backwards && t + dt <= time ) || (!backwards && t + dt >= time)
            dt = time - t 
            is_done = true
        end

        if !reuse_acceleration
            acc = f_acc(pos, vel, t)
        end

        pos_new, vel_new = step_kdk(pos, vel, acc, f_acc, dt, t)

        push!(positions, pos_new)
        push!(velocities, vel_new)
        push!(accelerations, acc)
        t = times[i] + dt
        push!(times, t)

        if is_done
            break
        end
    end

    positions_matrix = hcat(positions...)
    velocities_matrix = hcat(velocities...)
    accelerations_matrix = hcat(accelerations...)
    return Orbit(times=times, positions=positions_matrix, velocities=velocities_matrix, acceleration=accelerations_matrix)
end



function step_kdk(position::Vector{Float64}, velocity::Vector{Float64}, acceleration, f, dt::Real, t=0)
    vel_h = vel + dt/2 * acceleration
    pos_new = position + dt * vel_h
    acc = f(position, velocity, t + dt)
    vel_new = vel_h + 1/2 * dt * acc

    return pos_new, vel_new, acc
end


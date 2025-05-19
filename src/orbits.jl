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

    "Tidal stress at each time. While stress is technically a 3x3 matrix, here, we store the components as xx, yy, zz, xy, yz, zx"
    stresses::Union{Matrix{F}, Nothing} = nothing

    "pericentre with respect to origin"
    pericenter::F = minimum(calc_r(positions))

    "apocentre with respect to origin"
    apocenter::F = maximum(calc_r(positions))
end

positions(a::Orbit) = a.positions
velocities(a::Orbit) = a.velocities
accelerations(a::Orbit) = a.accelerations
stresses(a::Orbit) = a.stresses
pericenter(a::Orbit) = a.pericenter
apocenter(a::Orbit) = a.apocenter


Base.length(a::Orbit) = length(a.times)

""" Loads an orbit from a CSV file. Expected columns are `time`, `x`, `y`, `z`, `v_x`, `v_y`, `v_z`. """
function Orbit(filename::String)
    df = CSV.read(filename)
    positions = hcat(df.x, df.y, df.z)
    velocities = hcat(df.v_x, df.v_y, df.v_z)
    accelerations = haskey(df, "a_x") ? hcat(df.a_x, df.a_y, df.a_z) : nothing
    stresses = haskey(df, "da_xx") ? hcat(df.da_xx, df.da_yy, df.da_zz, df.da_xy, df.da_yz, df.da_zx) : nothing

    return Orbit(
        times = df.times,
        positions = positions,
        velocities = velocities,
        accelerations = accelerations,
        stresses = stresses,
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

    if a.stresses !== nothing
        df = hcat(df, DataFrame(
            da_xx = a.stresses[1, :],
            da_yy = a.stresses[2, :],
            da_zz = a.stresses[3, :],
            da_xy = a.stresses[4, :],
            da_yz = a.stresses[5, :],
            da_zx = a.stresses[6, :]
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
    @assert isapprox(a.time, b.time, rtol=1e-6)

    return Orbit(times=a.times, positions=a.positions - b.positions, velocities=a.velocities - b.velocities)
end


"""
    reverse(a::Orbit)

Reverses the order in time of an orbit
"""
function reverse(a::Orbit)
    time = reverse(a.time)
    positions = reverse(a.positions, dims=2)
    velocities = reverse(a.velocities, dims=2)
    acceleration = a.accelerations === nothing ? nothing : reverse(a.accelerations, dims=2)
    stresses = a.stresses === nothing ? nothing : reverse(a.stresses, dims=2)
    return Orbit(times=time, positions=positions, velocities=velocities, accelerations=accelerations, stresses=stresses)
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
- `acceleration::Function`: the acceleration function
- `dt_max::Real=0.1`: the maximum timestep
- `dt_min::Real=0.001`: the minimum timestep, will exit if timestep is below this
- `time::Real=-10/T2GYR`: the time to integrate to
- `timebegin::Real=0`: the time to start integrating from
- `timestep::Symbol=:adaptive`: the timestep to use, either `:adaptive` or a real number
- `η::Real=0.01`: the adaptive timestep parameter

"""
function leap_frog(gc, acceleration; 
        dt_max=0.1, dt_min=0.001, timebegin=0, time=-10/T2GYR, 
        timestep=:adaptive, η=0.01
    )

    if timestep isa Real
        dt_min = timestep
    end

    t = timebegin

    # setup matricies

    Nt = round(Int, abs(time / dt_min))
    positions = Vector{Vector{Float64}}()
    velocities = Vector{Vector{Float64}}()
    accelerations = Vector{Vector{Float64}}()
    times = Float64[]

    push!(positions, [gc.x, gc.y, gc.z])
    push!(velocities, [gc.v_x, gc.v_y, gc.v_z] / V2KMS)
    push!(accelerations, acceleration(positions[1], velocities[1], t))
    push!(times, 0.)
    is_done = false
    backwards = time < 0

    for i in 1:Nt
        pos = positions[i]
        vel = velocities[i]
        acc = acceleration(pos, vel, t)
        
        if timestep == :adaptive
            dt = min(sqrt(η / calc_r(acc)), dt_max)
        elseif timestep isa Real
            dt = timestep
        end
        
        if backwards
            dt *= -1
        end
        if abs(dt) < dt_min
            @warn "timestep below minimum timestep"
            break
        end

        if (backwards && t + dt <= time ) || (!backwards && t + dt >= time)
            dt = time - t
            is_done = true
        end

        pos_new, vel_new = step_dkd(pos, vel, acc, dt, t)

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



function step_kdk(position::Vector{Float64}, velocity::Vector{Float64}, acceleration, dt::Real, t=0)
    acc = acceleration(pos, vel, t)
    vel_h = vel + 1/2 * dt * acc
    pos_new = pos + dt*vel_h
    acc = acceleration(pos_new, vel_h, t + dt/2)
    vel_new = vel_h + 1/2 * dt * acc

    return pos_new, vel_new
end


"""
    step_dkd(position::Vector{Float64}, velocity::Vector{Float64}, acceleration, dt::Real, t=0)

Computes a single step of the leap frog method (drift-kick-drift) for a single particle.
Acceleration should be a function of (position, velocity, time)
"""
function step_dkd(position::Vector{Float64}, velocity::Vector{Float64}, acceleration, dt::Real, t=0)
    pos_h = pos + 1/2*dt* vel
    acc = acceleration(pos_new, vel_h, t + dt/2)
    vel_new = vel + dt * acc

    pos_new = pos_h + 1/2*dt* vel_new

    return pos_new, vel_new
end

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

angular_momenta(a::Orbit) = angular_momenta(positions(a), velocities(a))

radii(a::Orbit) = radii(a.positions)
speeds(a::Orbit) = radii(a.velocities)

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
    getindex(a::Orbit, idx)

Return an orbit only with timesteps selected by the given range or filter
"""
function Base.getindex(a::Orbit, idx)
    time = a.times[idx]
    positions = a.positions[:, idx]
    velocities = a.velocities[:, idx]
    accelerations = a.accelerations === nothing ? nothing : getindex(a.accelerations, :, idx)
    return Orbit(times=time, positions=positions, velocities=velocities, accelerations=accelerations)
end

"""
    resample(a::Orbit, time::AbstractVector)
Resamples an orbit to a new time array
"""
function resample(a::Orbit, time::AbstractVector{<:Real})
    if a.times[2] < a.times[1]
        a = reverse(a)
    end
    x = LilGuys.lerp(a.times, a.positions[1, :]).(time)
    y = LilGuys.lerp(a.times, a.positions[2, :]).(time)
    z = LilGuys.lerp(a.times, a.positions[3, :]).(time)
    v_x = LilGuys.lerp(a.times, a.velocities[1, :]).(time)
    v_y = LilGuys.lerp(a.times, a.velocities[2, :]).(time)
    v_z = LilGuys.lerp(a.times, a.velocities[3, :]).(time)

    return Orbit(times=time, positions=[x y z]', velocities=[v_x v_y v_z]')
end


"""
    apfrog(acceleration, gc; params...)

Computes the orbit of a Galactocentric object using the leapfrog integration method.

Parameters:
- `gc::Galactocentric`: the initial conditions
- `acceleration::Function`: the acceleration function. should take a position vector, velocity vector, and the current time (for maximum flexibility).
- `dt_max::Real=0.1`: the maximum timestep
- `dt_min::Real=0.001`: the minimum timestep, will exit if timestep is below this
- `timerange::Real=(0, -10/T2GYR)`: the initial and time to integrate to
- `timestep::Symbol=:adaptive`: the timestep to use, either `:adaptive` or a real number
- `η::Real=0.003`: the adaptive timestep parameter. Adaptive timestep is based on sqrt(η*r/a) where a is the magnitude of the current acceleration and r is the distance from the origin.

"""
function leapfrog(f_acc, gc::Galactocentric; 
        dt_max=0.1, dt_min=0.001, timerange=(0, -10/T2GYR), 
        timestep=:adaptive, η=0.003, reuse_acceleration=true,
        method = step_kdk,
    )


    if timestep isa Real
        dt_min = timestep
    end

    time_i, time_f = timerange
    integration_time = time_f - time_i
    t = time_i

    Nt = round(Int, abs(integration_time / dt_min))

    # setup matricies
    positions = Vector{Vector{Float64}}()
    velocities = Vector{Vector{Float64}}()
    accelerations = Vector{Vector{Float64}}()
    times = Float64[]

    push!(positions, [gc.x, gc.y, gc.z])
    push!(velocities, [gc.v_x, gc.v_y, gc.v_z] / V2KMS)
    push!(accelerations, f_acc(positions[1], velocities[1], t))
    push!(times, t)
    is_done = false
    backwards = time_f < time_i

    for i in 1:Nt
        pos = positions[i]
        vel = velocities[i]
        acc = accelerations[i]
        
        if timestep == :adaptive
            dt = min(sqrt(η * radii(pos) / radii(acc)), dt_max)
        elseif timestep isa Real
            dt = timestep
        end
        
        if backwards
            dt *= -1
        end
        if abs(dt) < dt_min
            @warn "timestep below minimum timestep"
            is_done = true
            dt = dt_min
        end

        # last step ends at end of range
        if (backwards && t + dt <= time_f ) || (!backwards && t + dt >= time_f)
            dt = time_f - t 
            is_done = true
        end

        if !reuse_acceleration
            acc = f_acc(pos, vel, t)
        end

        pos_new, vel_new, acc = method(pos, vel, acc, f_acc, dt, t)

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
    return Orbit(times=times, positions=positions_matrix, velocities=velocities_matrix, accelerations=accelerations_matrix)
end




function step_kdk(position::AbstractVector{<:Real}, velocity::AbstractVector{<:Real}, acceleration::AbstractVector{<:Real}, f::Function, dt::Real, t::Real=0)

    vel_h = velocity + acceleration * dt/2 # half kick
    pos_new = position + vel_h * dt # full drift
    acc_new = f(pos_new, velocity, t + dt) # acc_new (assuming vel. independ)
    vel_new = vel_h + acc_new * dt/2 # half kick

    return pos_new, vel_new, acc_new
end

function step_dkd(position::AbstractVector{<:Real}, velocity::AbstractVector{<:Real}, acceleration::AbstractVector{<:Real}, f::Function, dt::Real, t::Real=0)
    pos_h = position + velocity * dt/2  # half drift
    acc_h = f(pos_h, velocity, t + dt/2)
    vel_new = velocity + acc_h * dt # kick
    pos_new = pos_h +  vel_new * dt/2 # half drift

    return pos_new, vel_new, acc_h
end



"""
    peri_apo(orbit::Orbit)

Calculate the global peri and apocentre.
"""
function peri_apo(orbit)
    rs = radii(orbit)
    return minimum(rs), maximum(rs)
end



"""
    all_peris_apos(orbit)
    all_peris_apos(times, radii)

Compute all the pericentres and apocentres in the orbit

If an orbit is passed, assumes methods times(orbit) and radii(orbit) are defined.

Returns a 5-tuple of the pericentres, indices of pericentres, apocentres, indices of apocentres, and the pericentre error maximum.
"""
function all_peris_apos(ts::AbstractVector{<:Real}, rs::AbstractVector{<:Real})
    peris = Float64[]
    idx_peris = Int[]
    apos = Float64[]
    idx_apos = Int[]

    max_err = 0
    Nt = length(ts)

    for i in 2:Nt-1
        r = rs[i]
        last_r = rs[i-1]
        next_r = rs[i+1]

        if (last_r >= r) && (r <= next_r)
            # r is a local minimum
            push!(idx_peris, i)
            push!(peris, r)

            max_err = max(max_err, last_r - r, next_r - r)
        elseif (last_r <= r) && (r >= next_r)
            # is a maximum
            push!(idx_apos, i)
            push!(apos, r)
            max_err = max(max_err, r - last_r, r - next_r)
        end
    end

    return peris, idx_peris, apos, idx_apos, max_err
end

function all_peris_apos(orbit)
    return all_peris_apos(times(orbit), radii(orbit))
end


"""
    last_time_peri(ts, rs)

Compute the last pericentre and time of. Returns a tuple of the time of last pericentre and the pericentre, or NaNs if no solution is found. 
"""
function last_time_peri(ts::AbstractVector{<:Real}, rs::AbstractVector{<:Real})
    time = NaN
    peri = NaN

    Nt = length(ts)

    for i in 2:Nt-1
        r = rs[i]
        last_r = rs[i-1]
        next_r = rs[i+1]

        if (last_r >= r) && (r <= next_r)
            # r is a local minimum
            time = ts[i]
            peri = r
            break
        end
    end

    return time, peri
end


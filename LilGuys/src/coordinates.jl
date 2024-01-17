using Printf

Base.@kwdef struct PhasePoint
    x::F
    y::F
    z::F
    v_x::F
    v_y::F
    v_z::F
end

function Base.getproperty(pp::PhasePoint, sym::Symbol)
    if sym == :position
        return [pp.x, pp.y, pp.z]
    elseif sym == :velocity
        return [pp.v_x, pp.v_y, pp.v_z]
    elseif sym ∈ fieldnames(PhasePoint)
        return getfield(pp, sym)
    else
        error("$(sym) is not a valid field")
    end
end



Base.@kwdef struct Observation
    ra::F
    dec::F
    distance::F
    pm_ra::F
    pm_dec::F
    radial_velocity::F
end

function Observation(skycoord::PyObject)
    return Observation(
        skycoord.ra[],
        skycoord.dec[],
        skycoord.distance[],
        skycoord.pm_ra_cosdec[],
        skycoord.pm_dec[],
        skycoord.radial_velocity[],
        )
end

function PhasePoint(pos::Vector{F}, vel::Vector{F})
    if length(pos) != 3
        error("position must be a 3-vector")
    end
    if length(vel) != 3
        error("velocity must be a 3-vector")
    end
    return PhasePoint(pos..., vel...)
end

function PhasePoint(skycoord::PyObject)
    pos = [skycoord.x[], skycoord.y[], skycoord.z[]]
    vel = [skycoord.v_x[], skycoord.v_y[], skycoord.v_z[]]
    pos ./= R0
    vel ./= V0
    return PhasePoint(pos, vel)
end

function to_galcen(obs::Observation)
    skycoord = to_astropy(obs)
    transformed_coord = skycoord.transform_to(galcen_frame)
    return PhasePoint(transformed_coord)
end


function to_sky(phase::PhasePoint)

    skycoord = to_astropy(phase)
    transformed_coord = skycoord.transform_to(geocen_frame)

    return Observation(transformed_coord)
end


function to_astropy(obs::Observation)
    skycoord = astropy_coords.SkyCoord(
        ra = obs.ra * u.degree,
        dec = obs.dec * u.degree,
        distance = obs.distance * u.kpc,
        radial_velocity = obs.radial_velocity * u.km/u.s,
        pm_ra_cosdec = obs.pm_ra * u.mas/u.yr,
        pm_dec = obs.pm_dec * u.mas/u.yr
        )  
    return skycoord
end

function to_astropy(phase::PhasePoint)
    pos = phase.position * R0
    vel = phase.velocity * V0
    skycoord = astropy_coords.SkyCoord(
        x = pos[1] * R_U,
        y = pos[2] * R_U,
        z = pos[3] * R_U,
        v_x = vel[1] * V_U,
        v_y = vel[2] * V_U,
        v_z = vel[3] * V_U,
        frame = galcen_frame
        )

    return skycoord
end


"""
Returns a list of observations based on snapshot particles
"""
function to_sky(snap::Snapshot; invert_velocity::Bool=false, verbose::Bool=false)
    observations = Observation[]

    for i in 1:length(snap)
        if verbose
            print("converting $(i)/($(length(snap))\r")
        end
        pos = snap.positions[:, i]
        vel = snap.velocities[:, i]
        if invert_velocity
            vel *=-1
        end
        phase = PhasePoint(pos, vel)
        obs = to_sky(phase)
        push!(observations, obs)
    end
    return observations
end


function rand_coord(obs::Observation, err::Observation)
    return Observation(
        ra = obs.ra,
        dec = obs.dec,
        pm_ra = obs.pm_ra + randn() * err.pm_ra,
        pm_dec = obs.pm_dec + randn() * err.pm_dec,
        radial_velocity = obs.radial_velocity + randn() * err.radial_velocity,
        distance = obs.distance + randn() * err.distance,
       )
end

function rand_coords(obs::Observation, err::Observation, N::Int)
    return [rand_coord(obs, err) for _ in 1:N]
end


function Base.show(io::IO, pp::PhasePoint)
    x = pp.position * R0
    v = pp.velocity * V0
    @printf io "(%4.2f, %4.2f, %4.2f) kpc, " x...
    @printf io "(%4.2f, %4.2f, %4.2f) km/s" v...
    return io
end


function Base.show(io::IO, obs::Observation)
    @printf io "observation at (%4.2f, %4.2f)"  obs.ra obs.dec
    return io
end





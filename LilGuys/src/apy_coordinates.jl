module AstropyCoordinates


using PyCall
const astropy_coords = PyNULL()

import LilGuys as lguys

const u = PyNULL()
const galcen_frame = PyNULL()
const geocen_frame = PyNULL()


function __init__()
    copy!(u, pyimport("astropy.units"))
    copy!(astropy_coords, pyimport("astropy.coordinates"))

    astropy_coords.galactocentric_frame_defaults.set("v4.0")

    copy!(galcen_frame, astropy_coords.Galactocentric()) 
    copy!(geocen_frame, astropy_coords.ICRS())
end

function PhasePoint(skycoord::PyObject)
    pos = [skycoord.x[], skycoord.y[], skycoord.z[]]
    vel = [skycoord.v_x[], skycoord.v_y[], skycoord.v_z[]]
    pos ./= R0
    vel ./= V0
    return lguys.PhasePoint(pos, vel)
end

function Observation(skycoord::PyObject)
    return lguys.Observation(
        skycoord.ra[],
        skycoord.dec[],
        skycoord.distance[],
        skycoord.pm_ra_cosdec[],
        skycoord.pm_dec[],
        skycoord.radial_velocity[],
        )
end

function to_sky(phase::lguys.PhasePoint)

    skycoord = to_astropy(phase)
    transformed_coord = skycoord.transform_to(geocen_frame)

    return lguys.Observation(transformed_coord)
end


function to_astropy(obs::lguys.Observation)
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

function to_astropy(phase::lguys.PhasePoint)
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

function to_galcen(obs::lguys.Observation)
    skycoord = to_astropy(obs)
    transformed_coord = skycoord.transform_to(galcen_frame)
    return transformed_coord
end
end

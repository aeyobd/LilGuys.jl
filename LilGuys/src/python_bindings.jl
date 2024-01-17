using PyCall


const u = PyNULL()
const astropy_coords = PyNULL()
const galcen_frame = PyNULL()
const geocen_frame = PyNULL()
const galpy = PyNULL()
const galpy_covert = PyNULL()

function __init__()
    copy!(u, pyimport("astropy.units"))
    copy!(astropy_coords, pyimport("astropy.coordinates"))

    astropy_coords.galactocentric_frame_defaults.set("v4.0")

    copy!(galcen_frame, astropy_coords.Galactocentric()) 
    copy!(geocen_frame, astropy_coords.ICRS())

    pyimport("galpy.potential")
    pyimport("galpy.orbit")
    copy!(galpy, pyimport("galpy"))
    copy!(galpy_covert, pyimport("galpy.util.conversion"))

    copy!(M_U, u.Msun)
    copy!(T_U, u.Gyr)
    copy!(R_U, u.kpc)
    copy!(V_U, u.km/u.s)


    copy!(M_AP, M_U * M0)
    copy!(T_AP, T_U * T0)
    copy!(R_AP, R_U * R0)
    copy!(V_AP, V_U * V0)
end


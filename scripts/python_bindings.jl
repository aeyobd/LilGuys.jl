const galpy = PyNULL()
const galpy_covert = PyNULL()

function __init__()

    pyimport("galpy.potential")
    pyimport("galpy.orbit")
    copy!(galpy, pyimport("galpy"))
    copy!(galpy_covert, pyimport("galpy.util.conversion"))

    copy!(M_U, u.Msun)
    copy!(T_U, u.Gyr)
    copy!(R_U, u.kpc)
    copy!(V_U, u.km/u.s)
end


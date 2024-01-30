import LilGuys as lguys
using PyCall

import LinearAlgebra: transpose!

galpy = pyimport("galpy")
u = pyimport("astropy.units")
pyimport("galpy.potential")
pyimport("galpy.orbit")

F = lguys.F
R_U = u.kpc
V_U = u.km / u.s
M_U = u.M_sun
M0 = lguys.M0
R0 = lguys.R0
V0 = lguys.V0

@kwdef struct MWParams
    M_thin_disk::F = 5.9
    a_thin_disk::F = 3.9
    b_thin_disk::F = 0.31
    M_thick_disk::F = 2
    a_thick_disk::F = 4.4
    b_thick_disk::F = 0.92
    M_bulge::F = 2.1
    a_bulge::F = 1.3
    M_halo::F = 79.28
    a_halo::F = 20.2
end

function MWPotential(params::MWParams = MWParams())
    pot = galpy.potential

    Φ_bulge = pot.HernquistPotential(
        amp=params.M_bulge * M_U * M0,
        a=params.a_bulge * R_U * R0)
	Φ_halo = pot.NFWPotential(
        amp=params.M_halo * M0 * M_U,
        a=params.a_halo * R0 * R_U)
	Φ_thick_disk = pot.MiyamotoNagaiPotential(
        amp=params.M_thick_disk * M0 * M_U,
        a=params.a_thick_disk * R0 * R_U,
        b=params.b_thick_disk * R0 * R_U)
	Φ_thin_disk = pot.MiyamotoNagaiPotential(
        amp=params.M_thin_disk * M0 * M_U,
        a=params.a_thin_disk * R0 * R_U, 
        b=params.b_thin_disk * R0 * R_U)


	Φ_mw = Φ_bulge + Φ_halo + Φ_thick_disk + Φ_thin_disk
    return Φ_mw
end

function eval_Φ(Φ::PyObject, x, y, z)
    return galpy.potential.evaluatePotentials(Φ, x * R0 * R_U, y * R0 * R_U, z * R0 * R_U)

end

function to_cylendrical(phase::lguys.Galactocentric)
    # galactocentric is in natural units
    x, y, z = phase.position
    v_x, v_y, v_z = phase.velocity

    R = sqrt(x^2 + y^2) 

    vR = (x*v_x + y*v_y) / R * V_U 
    vT = (x*v_y - y*v_x) /R * V_U
    vz = phase.v_z * V_U

    phi = atan(phase.y, phase.x) * u.rad

    R *= R_U
    z *= R_U
    return R, vR, vT, z, vz, phi
end

function calc_orbit(phase::lguys.Galactocentric; 
        times=LinRange(0, 10, 10000), Φ=MWPotential())
    o = galpy.orbit.Orbit(to_cylendrical(phase))

    t = times * u.Gyr

    o.integrate(t, Φ, method="odeint")
    positions = hcat(o.x(t), o.y(t), o.z(t))
    velocities = hcat(o.vx(t), o.vy(t), o.vz(t))

    positions = Matrix(transpose(positions))
    velocities = Matrix(transpose(velocities))
    Φs = eval_Φ(Φ, o.x(t), o.y(t), o.z(t))

    return t, positions, velocities, Φs
end


function V_circ(Φ::PyObject, r)
    return galpy.potential.vcirc(Φ, r * R0 * u.kpc) 
end


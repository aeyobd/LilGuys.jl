
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

function calc_orbit(phase::PhasePoint; 
        times=LinRange(0, 10, 10000), Φ=MWPotential())
    R = sqrt(phase.x^2 + phase.y^2) * R_U * R0
    z = phase.z * R_U * R0
    vR = phase.v_x * V_U * V0
    vT = phase.v_y * V_U * V0
    vz = phase.v_z * V_U * V0
    phi = atan(phase.y / phase.x) * u.rad

    o = galpy.orbit.Orbit([R,vR,vT,z,vz,phi])

    t = times * u.Gyr

    o.integrate(t, Φ, method="odeint")
    x = o.x(t)
    y = o.y(t)
    z = o.z(t)
    vx = o.vx(t)
    vy = o.vy(t)
    vz = o.vz(t)

    Φs = eval_Φ(Φ, x, y, z)
    return t, x, y, z, vx, vy, vz, Φs
end


function V_circ(Φ::PyObject, r)
    return galpy.potential.vcirc(Φ, r * R0 * u.kpc) 
end

module AgamaExt


import Agama
using LilGuys

import LilGuys: AgamaPotential, sample_potential, agama_orbit, leapfrog



function AgamaPotential(nfw::LilGuys.TruncNFW)
    r_s = nfw.r_s
    M_s = nfw.M_s
    rho_s = M_s / (4π * r_s^3)
    r_t = nfw.r_t

    return Agama.Potential(type="Spheroid", densityNorm=rho_s,
                gamma = 1, beta=3, alpha=1, scaleRadius=r_s, 
            outerCutoffRadius = r_t, cutoffStrength=1)
end



function sample_potential(nfw::LilGuys.TruncNFW, N::Integer=10000; kwargs...)
    pot = AgamaPotential(nfw)
    df = Agama.DistributionFunction(;type="QuasiSpherical", potential=pot._py, kwargs...)
    gm = Agama.GalaxyModel(pot, df)
    pos, vel, mass = Agama.sample(gm, N)

    return LilGuys.Snapshot(pos, vel, mass)
end



function LilGuys.Orbit(o::Agama.Orbit)
    return LilGuys.Orbit(times=Agama.times(o), positions=Agama.positions(o), velocities=Agama.velocities(o))
end


"""
    agama_orbit(pot::Agama.Potential, ic; kwargs...)

Compute the orbit given an Agama potential and a LilGuys.CoordinateFrame or 
vector of coordinates.
Returns an orbit or a vector of orbits respectively.
"""
function agama_orbit(pot::Agama.Potential, coords::AbstractVector{<:LilGuys.CoordinateFrame}; agama_units=Agama.current_units(), kwargs...)
    coords_i = LilGuys.transform.(LilGuys.Galactocentric, coords)
    pos_i = hcat(LilGuys.position.(coords_i)...)
    vel_i = hcat(LilGuys.velocity.(coords_i)...) ./ V2KMS

    o = Agama.orbit(pot, pos_i, vel_i, agama_units; kwargs...)

    return LilGuys.Orbit.(o)
end

function agama_orbit(pot::Agama.Potential, coord::LilGuys.CoordinateFrame; agama_units=Agama.current_units(), kwargs...)
    coords_i = LilGuys.transform(LilGuys.Galactocentric, coord)
    pos_i = LilGuys.position(coords_i)
    vel_i = LilGuys.velocity(coords_i) ./ V2KMS

    o = Agama.orbit(pot, pos_i, vel_i, agama_units; kwargs...)

    return LilGuys.Orbit(o)
end


"""
   leapfrog(pot::Agama.Potential, coord; kwargs...)

    Integrate in given potential. See alternate docstring for `leapfrog`, only changed option is `agama_units` is passed to `Agama`
"""
function leapfrog(pot::Agama.Potential, coord::LilGuys.Galactocentric; agama_units=Agama.current_units(), kwargs...)
    f_acc = (pos, vel, t) -> Agama.acceleration(pot, pos, agama_units, t=t)
    return LilGuys.leapfrog(f_acc, coord; kwargs...)
end



end # module

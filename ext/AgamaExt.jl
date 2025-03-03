module AgamaExt


import Agama
using LilGuys

import LilGuys: AgamaPotential, sample_potential



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



end # module

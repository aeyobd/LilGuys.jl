module LilGuys

using Requires


export Snapshot
export Output
export ICRS, PhasePoint
export to_galcen, to_sky
export save

export calc_r, calc_v


# profile tools
export AbstractProfile, NFW
export calc_ρ, calc_M, calc_r_circ_max, calc_v_circ_max, calc_v_circ

export M2MSUN, R2KPC, V2KMS, T2GYR
export ICRS, HelioRest, Galactocentric, transform



include("units.jl")
include("interface.jl")
using .Interface


include("utils.jl")

include("coordinates.jl")

include("io.jl")
include("snapshot.jl")
include("output.jl")

include("spherical.jl")
include("coord_trans.jl")   

include("physics.jl")
include("gravity.jl")

include("analytic_profiles.jl")
include("nfw.jl")
include("scaling_relations.jl")

include("project.jl")
include("density_3d.jl")
include("density_2d.jl")
include("density_3d_star.jl")

include("centres/Centres.jl")



function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("plots.jl")
end


using .Centres

end # module LilGuys

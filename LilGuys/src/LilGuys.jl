module LilGuys

export Snapshot
export Output
export Observation, PhasePoint
export to_galcen, to_sky
export save

export calc_r, calc_v

include("units.jl")
include("utils.jl")

include("coordinates.jl")

include("snapshot.jl")
include("output.jl")

include("coord_trans.jl")   

include("physics.jl")
include("gravity.jl")
include("analytic_profiles.jl")
include("density_3d.jl")
include("density_utils.jl")

include("centres/Centres.jl")


using Requires
function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("plots.jl")
end

using .Centres

end # module LilGuys

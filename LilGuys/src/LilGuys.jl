module LilGuys

export Snapshot
export Output
export Observation, PhasePoint
export to_galcen, to_sky
export save


include("units.jl")
include("utils.jl")

include("coordinates.jl")

include("io/snapshot.jl")
include("io/output.jl")

include("coord_trans.jl")   

include("physics.jl")
include("gravity.jl")
include("profile.jl")


include("centres/static_centres.jl")
include("centres/centre_fuzzy.jl")
include("centres/shrinking_spheres.jl")
include("centres/centre_output.jl")


using Requires
function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plots.jl")
end

end # module LilGuys

module LilGuys

export Snapshot
export Output
export Observation, PhasePoint
export to_galcen, to_sky
export save


include("units.jl")
include("utils.jl")

include("coordinates.jl")

include("snapshot.jl")
include("output.jl")

include("coord_trans.jl")   

include("physics.jl")
include("gravity.jl")
include("profile.jl")

include("centres/Centres.jl")


using Requires
function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plots.jl")
end

using .Centres

end # module LilGuys

module LilGuys

export Snapshot
export Output
export Observation, PhasePoint
export to_galcen, to_sky
export save


include("units.jl")
include("utils.jl")
# include("particle.jl")
include("coordinates.jl")

include("snapshot.jl")
include("phys_quantities.jl")
include("coord_trans.jl")   
include("gravity.jl")
include("output.jl")
include("profile.jl")

include("centre_fuzzy.jl")
include("centre_shrinking_spheres.jl")
include("centre_density.jl")

using Requires
function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plots.jl")
end

end # module LilGuys

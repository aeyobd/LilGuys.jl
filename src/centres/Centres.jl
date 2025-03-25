module Centres
    
    export Centre, SS_State, StaticState
    export calc_centre, calc_centres, calc_centre!
    export shrinking_spheres

    using ..LilGuys
    import ..centroid, ..centroid_err, ..radii, ..F, ..Snapshot, ..Output, ..quantile, ..mean
    import ..potential_spherical_discrete, ..potential_spherical_func
    import ..specific_energy
    import ..OptVector

    include("static_centres.jl")
    include("centre_output.jl")
    include("shrinking_spheres.jl")
    include("centre_fuzzy.jl")

end

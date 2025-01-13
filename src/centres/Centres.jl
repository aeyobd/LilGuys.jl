module Centres
    
    export SS_State, calc_centre, calc_centres, StaticState
    export calc_centre!
    export shrinking_spheres

    using ..LilGuys
    import ..centroid, ..centroid_err, ..calc_r, ..F, ..Snapshot, ..Output, ..quantile, ..mean
    import ..calc_radial_discrete_Φ
    import ..calc_E_spec
    import ..OptVector

    include("static_centres.jl")
    include("centre_output.jl")
    include("shrinking_spheres.jl")
    #include("centre_fuzzy.jl")

end

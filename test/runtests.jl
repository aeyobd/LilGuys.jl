include("setup.jl")


tests = ["units", "utils", 
         "measurements",
         "io", "interface",
         "snapshot", "output",
         "spherical", "coordinates", "coord_trans",
         "surface_density", "density_profile", "mass_profile",
         "stellar_density_3d", "project",
         "nfw",
         "physics", "gravity", 
         "scaling_relations",
         "centre_static", "analytic_profiles", 
         "centre_output",
         "shrinking_spheres", 
         "plots",
        ]

for test in tests
    @testset "$test" begin
        include("$(test)_tests.jl")
    end
end

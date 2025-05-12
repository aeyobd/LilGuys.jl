include("setup.jl")


tests = ["units", "utils", 
         "measurements",
         "io", "interface",
         "snapshot", "output",
         "spherical", "coordinates", "coord_trans",
         "density_2d", "density_3d", "mass_profile_3d"
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

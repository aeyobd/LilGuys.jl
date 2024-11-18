include("setup.jl")


tests = ["units", "utils", 
         "io", "interface",
         "snapshot", "output",
         "spherical", "coordinates", "coord_trans",
         "density_2d", "density_3d", "density_3d_star",
         "nfw",
         "physics", "gravity", 
         "centre_static", "analytic_profiles", 
         "shrinking_spheres",
        ]

for test in tests
    @testset "$test" begin
        include("$(test)_tests.jl")
    end
end

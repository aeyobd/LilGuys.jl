include("setup.jl")


tests = ["units", "utils", 
         "io",
         "snapshot", "profile", 
         "spherical", "coordinates", "coord_trans",
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

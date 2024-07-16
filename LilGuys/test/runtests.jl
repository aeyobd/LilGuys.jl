include("setup.jl")


tests = ["units", "utils", 
         "io",
         "snapshot", "profile", 
         "coordinates", "coord_trans",
         "physics", "gravity", 
         "centre_static", "profiles", "spherical",
         "shrinking_spheres",
        ]

for test in tests
    @testset "$test" begin
        include("$(test)_tests.jl")
    end
end

include("setup.jl")


tests = ["utils", "snapshot", "profile", "coordinates", "phys_quantities", "gravity", "centre_shrinking_spheres"]

for test in tests
    @testset "$test" begin
        include("$(test)_tests.jl")
    end
end

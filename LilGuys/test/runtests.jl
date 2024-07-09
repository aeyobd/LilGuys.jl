include("setup.jl")


tests = ["utils", "snapshot", "profile", "coordinates", "physics", "gravity", "centre_static", "profiles", "spherical"]

for test in tests
    @testset "$test" begin
        include("$(test)_tests.jl")
    end
end

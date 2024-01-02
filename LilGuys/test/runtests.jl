include("setup.jl")


tests = ["utils", "hdf5_utils", "snapshot", "fuzzy_snapshot", "profile", "coordinates", "phys_quantities", "gravity", "centre"]

for test in tests
    @testset "$test" begin
        include("$(test)_tests.jl")
    end
end

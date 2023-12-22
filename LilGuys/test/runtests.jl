using Test
using LilGuys

@testset "hdf5 utils" begin
    include("hdf5_utils_tests.jl")
end

@testset "snapshot " begin
    include("snapshot_tests.jl")
end

@testset "profile" begin
    include("profile_tests.jl")
end

@testset "coordinates" begin
    include("coord_tests.jl")
end


tdir = mktempdir()

function struct_approx_equal(a, b, rtol=1e-8, atol=1e-8)
    for field in fieldnames(typeof(a))
        v1 = getfield(a, field)
        v2 = getfield(b, field)

        if (v1 isa AbstractFloat) || (v1 isa Array && eltype(v1) <: AbstractFloat)
            @test isapprox(v1, v2, rtol=rtol, atol=atol)
        else
            @test v1 == v2
        end
    end
end


@testset "write_hdf5_table & read_hdf5_table" begin
    @testset "simple" begin
        df = lguys.DataFrame(
            a = [1, 2, 8],
            funky_key = [4.0, -π, NaN],
            wow = ["a", "b", "c"],
        )
        df[!, Symbol("nasty test")] = [true, false, true]
        
        filename = joinpath(tdir, "test.hdf5")

        # write to fits
        lguys.write_hdf5_table(filename, df, overwrite=true)

        # read from fits
        df2 = lguys.read_hdf5_table(filename)

        # check if the two dataframes are equal
        @test size(df) == size(df2)
        @test Set(names(df2)) == Set(names(df))
        @test df.a == df2.a
    end


    @testset "exceptions" begin 
    end
end



@testset "write and read structs_from_hdf5" begin

    @testset "simple" begin
        @Base.kwdef struct TestStruct
            a::Int
            b::Float64
            c::String
            v::Vector{Float64}
        end

        s1 = TestStruct(1, 2.0, "a", [1.0, 2.0, 3.0])
        s2 = TestStruct(0, π, "b", [4.0, 5.0, -2.2])

        structs = ["1"=>s1, "2"=>s2]

        filename = joinpath(tdir, "test.hdf5")
        lguys.write_structs_to_hdf5(filename, structs)

        structs2 = lguys.read_structs_from_hdf5(filename, TestStruct)

        for i in 1:length(structs)
            struct_approx_equal(structs[i].second, structs2[i].second)
            @test structs[i].first == structs2[i].first
        end

    end

end

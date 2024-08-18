tdir = mktempdir()

@testset "write_fits & load_fits" begin

    @testset "simple" begin
        df = lguys.DataFrame(
            a = [1, 2, 8],
            funky_key = [4.0, -π, NaN],
            wow = ["a", "b", "c"],
        )
        df[!, Symbol("nasty test")] = [true, false, true]
        
        # write to fits
        lguys.write_fits(joinpath(tdir, "test.fits"), df)

        # read from fits
        df2 = lguys.load_fits(joinpath(tdir, "test.fits"))


        # check if the two dataframes are equal
        @test size(df) == size(df2)
        @test Set(names(df2)) == Set(names(df))
        @test df.a == df2.a
    end

    @testset "exceptions" begin 
        df = lguys.DataFrame(
            α = [1,2]
        )

        @test_throws ArgumentError lguys.write_fits(joinpath(tdir, "test.fits"), df)
    end


end


@testset "write_hdf5_table & load_hdf5_table" begin
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
        df2 = lguys.load_hdf5_table(filename)

        # check if the two dataframes are equal
        @test size(df) == size(df2)
        @test Set(names(df2)) == Set(names(df))
        @test df.a == df2.a
    end


    @testset "exceptions" begin 
    end
end

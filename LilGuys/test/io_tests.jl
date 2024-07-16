@testset "write and retrieve fits" begin
    # create sample dataframe
    # DataFrame is from DataFrames.jl

    df = lguys.DataFrame(
        a = [1, 2, 8],
        funky_key = [4.0, -π, NaN],
        wow = ["a", "b", "c"],
    )
    df[!, Symbol("nasty test")] = [true, false, true]
    
    # write to fits
    lguys.write_fits("test.fits", df)

    # read from fits
    df2 = lguys.load_fits("test.fits")


    # check if the two dataframes are equal
    @test size(df) == size(df2)
    @test Set(names(df2)) == Set(names(df))
    @test df.a == df2.a

end



# these tests are written mostly for basic functionality
# and to ensure the interface is as documented.
# Correctness is assumed since these are functions 
# from other packages (

@testset "stats" begin
    x = [1, -0.1, 0.3]
    w = [0.6, 0, 0.3]

    @test lguys.mean(x) ≈ 0.4
    @test lguys.mean(x, w) ≈ 0.69 / 0.9

    @test lguys.std(x) ≈ sqrt(((1-0.4)^2 + (-0.1-0.4)^2 + (0.3-0.4)^2)/2)
    @test lguys.std(x, w) ≈ sqrt( ((1-0.7666666666666666)^2 * 0.6 +(0.3 - 0.7666666666666666)^2 * 0.3) / 0.9)


    # quantile
    #
    # weighted quantile
end


@testset "histogram" begin


end

@testset "bins_default" begin


end

@testset "effective_sample_size" begin

end



@testset "integrate" begin

end


@testset "filter_nans" begin

end



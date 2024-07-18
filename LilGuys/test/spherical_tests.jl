@testset "angular_distance simple cases" begin
    @test lguys.angular_distance(0, 0, 0, 0) ≈ 0
    @test lguys.angular_distance(0, 0, 0, 90) ≈ 90
    @test lguys.angular_distance(0, 0, 90, 0) ≈ 90
    @test lguys.angular_distance(0, 0, -90, 0) ≈ 90
    @test lguys.angular_distance(0, 0, 0, -90) ≈ 90

end


@testset "angular_distance out of bound" begin


end


@testset "angular_distance identity" begin
    # tricy almost same test
    @test lguys.angular_distance.([14.479505360741172], 
                                  [-33.99473291943769], [14.479505360741172], [-33.99473291943769]) == [0.]



    N = 1000
    ra = rand(N) * 360
    dec = rand(N) * 180 .- 90

    @test lguys.angular_distance.(ra, dec, ra, dec) ≈ zeros(N)  atol=1e-4
end




@testset "to_tangent" begin
    @test false broken = true
end


@testset "unit_vector" begin
    @test false broken = true
end


@testset "caartesian_to_sky" begin
    @test false broken = true
end


@testset "Rx" begin
    @test false broken = true
end

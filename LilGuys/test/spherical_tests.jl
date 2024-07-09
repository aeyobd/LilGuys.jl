@testset "angular_distance" begin
    @test lguys.angular_distance(0, 0, 0, 0) ≈ 0
    @test lguys.angular_distance(0, 0, 0, 90) ≈ 90
    @test lguys.angular_distance(0, 0, 90, 0) ≈ 90
    @test lguys.angular_distance(0, 0, -90, 0) ≈ 90
    @test lguys.angular_distance(0, 0, 0, -90) ≈ 90

    # tricy almost same test
    @test lguys.angular_distance.([14.479505360741172], 
                                  [-33.99473291943769], [14.479505360741172], [-33.99473291943769]) == [0.]
end

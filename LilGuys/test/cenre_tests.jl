@testset "centroid" begin
    x = [1 2 3; 5 6 7]

    expected = [3 4 5;]

    cen = centroid(x)
    actual = [a for a in cen.pos]

    @test actual ≈ expected
end

@testset "centroid weights" begin
    x = [1 2 3; 4 5 6]
    w = [1 ; 2]

    expected = [3 4 5;]

    cen = centroid(x, w)
    actual = [a for a in cen.pos]

    @test actual ≈ expected

end



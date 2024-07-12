
@testset "centre of mass" begin

end

@testset "potential weighted" begin
end


@testset "centroid snap" begin
    N = 10
    pos = transpose(hcat(ones(N), zeros(N), -ones(N)))
    vel = transpose(hcat(zeros(N), 2*ones(N), zeros(N)))
    m = lguys.ConstVector(1., N)
    snap = lguys.Snapshot(pos, vel, m)

    # cen = lguys.centroid(snap)

    @test cen.pos ≈ [1., 0., -1.] skip=true
    @test cen.vel ≈ [0., 2., 0.] skip=true
end



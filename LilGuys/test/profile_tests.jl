
@testset "sort_by_r" begin
    N = 100
    m = randn(N)
    r = 10 .^ randn(N)
    cen = [0.1, 0.2, 0.4]

    pos = r' .* lguys.rand_unit(N)  .+ cen
    vel = zeros(3, N)

    snap = lguys.Snapshot(pos, vel, m)
    snap.x_cen = cen

    sorted = lguys.sort_by_r(snap)

    r_end = lguys.calc_r(sorted)

    @test issorted(r_end)
    @test m[sortperm(r)] == sorted.masses

end

@testset "calc_ρ_hist" begin

end

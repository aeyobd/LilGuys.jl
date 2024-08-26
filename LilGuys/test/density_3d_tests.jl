



@testset "calc_ρ_hist" begin

end

@testset "calc_v_rad" begin
    @testset "basic" begin

        # simple directions
        pos = [ [-1,0,0]     [0, 2, 0]     [1/2, 3/2, -1] [π, √π, π^π]]
        vel = [ [π, 2, -2]  [0, 2.1, 0]   [0, 0, 0]      [0,0,0]]

        @test lguys.calc_v_rad(pos, vel) ≈ [-π, 2.1, 0, 0]

        # easy angles
        pos = [ [3, 4, 0] [0, 12/3, -5/3]]
        vel = [ [2, 0, 0] [1, 0, 0.1]]

        @test lguys.calc_v_rad(pos, vel) ≈ [2*3/5, 0.1 * -5/13]
    end


    @testset "random" begin
        N = 100
        pos = randn(3, N)

        # random orthoganal
        y = randn(3, N)
        vel = rand(N)' .* hcat([y[:, i] × pos[:, i] for i in 1:N]...)

        @test lguys.calc_v_rad(pos, vel) ≈ zeros(N) atol=1e-13

        # random parallel
        v = 10 .^ randn(N)
        vel = v' .* pos ./ lguys.calc_r(pos)'

        @test lguys.calc_v_rad(pos, vel) ≈ v

        pos = pos .* (10 .^ randn(N))'
        @test lguys.calc_v_rad(pos, vel) ≈ v

        # antiparallel
        vel = -vel
        @test lguys.calc_v_rad(pos, vel) ≈ -v
    end
end


@testset "v_circ_max_model" begin
    param = (1, 1)
    halo = lguys.NFW(v_circ_max=1.0, r_circ_max=1.0)

    @test lguys._v_circ_max_model(1, param) ≈ 1.0
    @test lguys._v_circ_max_model(1.01, param) < 1
    @test lguys._v_circ_max_model(0.99, param) < 1
    @test lguys._v_circ_max_model(-1, param) === NaN

    param = (2.23, √2)
    halo = lguys.NFW(v_circ_max=param[2], r_circ_max=param[1])
    @test lguys._v_circ_max_model(param[1], param) ≈ param[2]

    x = LinRange(0.1, 10, 100)
    y = lguys._v_circ_max_model(x, param)
    y2 = lguys.calc_v_circ.(halo, x)

    @test y ≈ y2

    @test lguys._v_circ_max_model(2.23, (NaN, NaN)) === NaN
    @test lguys._v_circ_max_model(2.23, (-1, 2.4)) === NaN
    @test lguys._v_circ_max_model(π, (π, -2.4)) ≈ -2.4
end


@testset "fit_v_r_circ_max" begin
    @testset "simple cases" begin
        halo = lguys.NFW(v_circ_max=1.0, r_circ_max=1.0)
        r = 10 .^ LinRange(-1, 1, 100)
        v = lguys.calc_v_circ.(halo, r)

        fit = lguys.fit_v_r_circ_max(r, v)
        @test fit.r_circ_max ≈ 1.0
        @test fit.v_circ_max ≈ 1.0
        @test fit.converged


        halo = lguys.NFW(v_circ_max=π/2, r_circ_max=√π)
        r = 10 .^ LinRange(-1, 1, 100)
        v = lguys.calc_v_circ.(halo, r)

        fit = lguys.fit_v_r_circ_max(r, v)
        @test fit.r_circ_max ≈ √π
        @test fit.v_circ_max ≈ π/2
        @test fit.converged
    end

    @testset "cannot converge" begin
        r = 10 .^ LinRange(-1, 1, 100)
        v = fill(2.23, 100)

        fit = lguys.fit_v_r_circ_max(r, v)
        @test !fit.converged
        @test fit.v_circ_max ≈ 2.23
    end

    @testset "exceptions" begin
        r = [1.0, 2.0, 3.0]
        v = [1.0, 2.0]

        @test_throws DimensionMismatch lguys.fit_v_r_circ_max(r, v)

        @test_throws ArgumentError lguys.fit_v_r_circ_max([1.], [1.])
    end
end



@testset "calc_v_circ" begin
end


@testset "calc_ρ_hist" begin

end

@testset "calc_M_in" begin
    @testset "simple_cases" begin
        N = 10
        snap = lguys.Snapshot(positions=zeros(3, N), velocities=zeros(3, N), masses=ones(N), index=1:N, header=lguys.make_default_header(1, N), Φs=-ones(N))

        @test lguys.calc_M_in(snap, 0.01) == N


        r = collect(1:10) .- 0.01
        snap = lguys.Snapshot(positions=r' .* lguys.rand_unit(N), velocities=zeros(3, N), masses=ones(N), index=1:N, header=lguys.make_default_header(1, N), Φs=-ones(N))

        expected = collect(1:10)
        actual = [lguys.calc_M_in(snap, i) for i in 1:N]
        @test expected ≈ actual
    end
end


@testset "calc_profile (integration)" begin


end


@testset "ObsProfile3D" begin

    @testset "constructor" begin

    end

    @testset "integration with NFW" begin

    end
end




@testset "to_sky" begin
    @testset "inverse" begin
        # can we transform coordinate back and recover original snapshot?
    end
end


@testset "to_gaia" begin
end


@testset "to_frame" begin
end

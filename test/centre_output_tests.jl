function make_out_test(x_expected, v_expected; N=100, r_s=0.3, v_s=0.05, dt=0.1)
    dir = mktempdir()
    mkdir(joinpath(dir, "out"))


    Nt = size(x_expected, 2)
    header = LilGuys.make_default_header(N, 0)
    times = LinRange(0, Nt*dt, Nt)

    for i in 1:Nt
        positions = LilGuys.rand_unit(N) .* r_s .+ x_expected[:, i]
        velocities = LilGuys.rand_unit(N) .* v_s .+ v_expected[:, i]
        masses = rand(N)
        index = collect(3i:N+3i - 1) .% N .+ 1 # shift indicies for consistency 
        snap = Snapshot(positions, velocities, masses,
                        index=index, header=header, time=times[i])

        snap.Φs = lguys.calc_radial_discrete_Φ(snap)
        LilGuys.save("$dir/out/snapshot_$i.hdf5", snap)
    end


    old_dir = pwd()
    cd(dir * "/out")

    h5open("combined.hdf5", "w") do file
        for i in 1:Nt
            HDF5.create_external(file, "/snap$(i-1)", "snapshot_$i.hdf5", "/")
        end
    end

    cd(old_dir)
    out = Output(dir)

    return out
end


@testset "calc centres" begin
    x_expected = [0.2  0.1
                 -0.9 0.5
                  1.5  1.1
                   ]
    v_expected = [0.6  0.1
                  3.2  0.5
                  -0.8 1.1
                   ]
    N = 100

    r_s =  0.3
    v_s = 0.05

    out = make_out_test(x_expected, v_expected, N=N, r_s=r_s, v_s=v_s)
    Nt = length(out)


    @testset "calc simple centres" begin
        state = LilGuys.Centres.StaticState

        cens = LilGuys.Centres.calc_centres(state, out)
        @test length(cens) == Nt 
        x_cen = hcat([c.position for c in cens]...)
        v_cen = hcat([c.velocity for c in cens]...)

        @test x_cen ≈ x_expected rtol=1e-1
        @test v_cen ≈ v_expected rtol=1e-2
    end

    @testset "calc ss centres" begin
        state = LilGuys.Centres.SS_State

        cens = LilGuys.Centres.calc_centres(state, out)
        @test length(cens) == Nt 
        x_cen = hcat([c.position for c in cens]...)
        v_cen = hcat([c.velocity for c in cens]...)

        @test x_cen ≈ x_expected rtol=1e-1
        @test v_cen ≈ v_expected rtol=1e-2
    end

end


let 
    # by intentially seting up the velocity in the oposite direction
    # the prior will cause most particles to be cut, so reinit state will matter
    x_expected = [0.2  0.1
                 -0.9 0.5
                  1.5  1.1
                   ]
    v_expected = [-0.3  0.1
                  -0.5  0.5
                  0.2  1.1
                   ]
    N = 100

    r_s =  0.3
    v_s = 0.05

    out = make_out_test(x_expected, v_expected, N=N, r_s=r_s, v_s=v_s, dt=1)
    Nt = length(out)

    @testset "reinit centre" begin
        state = LilGuys.Centres.SS_State

        cens = LilGuys.Centres.calc_centres(state, out, 
            reinit_state=false, r_max=0.3)

        cens_reinit = LilGuys.Centres.calc_centres(state, out,
            reinit_state=true, r_max=0.3)

        cen_f = LilGuys.Centres.calc_centre(LilGuys.Centres.SS_State, 
            out[end], r_max=0.3)

        @test cen_f.position ≈ cens_reinit[end].position 
        @test cen_f.velocity ≈ cens_reinit[end].velocity
        @test calc_r(cens[end].position, cens_reinit[end].position) > 1e-5
        @test calc_r(cens[end].velocity, cens_reinit[end].velocity) > 1e-5

    end
end

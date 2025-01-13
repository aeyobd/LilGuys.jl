let
    dir = mktempdir()
    mkdir(joinpath(dir, "out"))

    header = LilGuys.make_default_header(2, 0)

    x_expected = [0.2  0.1
                 -0.9 0.5
                  1.5  1.1
                   ]
    v_expected = [0.6  0.1
                  3.2  0.5
                  -0.8 1.1
                   ]
    N = 100
    Nt = size(x_expected, 2)

    r_s =  0.3
    v_s = 0.05


    for i in 1:Nt
        positions = LilGuys.rand_unit(N) .* r_s .+ x_expected[:, i]
        velocities = LilGuys.rand_unit(N) .* v_s .+ v_expected[:, i]
        masses = rand(N)
        index = collect(3i:N+3i - 1) .% N .+ 1 # shift indicies for consistency 
        snap = Snapshot(positions=positions, velocities=velocities, masses=masses, 
                        index=index, header=header)
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

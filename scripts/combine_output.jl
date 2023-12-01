using LilGuys

using Glob
using HDF5

function parse_args()
    if length(ARGS) != 2
        print("expected 2 arguments")
    end
    din, fout = ARGS
    return din, fout
end

function check_dir(fout)
    if ispath(fout)
        println("overwrite $(fout)? (y/N)")
        a = readline()
        if a != "y"
            error("file exists")
        end
        rm(fout)
    end
end

function get_filenames(din)
    filenames =  glob("snapshot*.hdf5", din)

    ms = [match(r"\d+", splitdir(f)[end]) for f in filenames]
    idx = [parse(Int, m.match) for m in ms]
    perm = sortperm(idx)

    return filenames[perm]
end

function new_file(func, filename, Nt, Np)
    h5open(filename, "w") do f
        create_group(f, "PartType1")

        vshape = (Nt, 3, Np)
        for name in ["Coordinates", "Velocities", "Acceleration"]
            create_dataset(f, "PartType1/" * name, F, vshape)
        end

        scalar_shape = (Nt, Np)
        for name in ["Potential", "ExtPotential"]
            create_dataset(f, "PartType1/" * name, F, scalar_shape)
        end

        create_dataset(f, "PartType1/ParticleIDs", Int, (Nt, Np))
        create_dataset(f, "Time", F, (Nt,))

        println("file opened")
        println()

        func(f)
    end
end


function combine_outputs(din, fout)
    filenames = get_filenames(din)
    snap1 = Snapshot(filenames[1])

    index1 = sort(snap1.index)
    Nt = length(filenames)
    Np = length(snap1)

    new_file(fout, Nt, Np) do f
        LilGuys.set_header!(f, snap1.header)
        for i in 1:Nt
            print("$(i)/$(Nt) \r")

            snap = Snapshot(filenames[i])
            perm = sortperm(snap.index)
            snap = snap[perm]
            if all(snap.index .!= index1)
                error("Mixmatched indicies")
            end
            for (sym, name) in LilGuys.h5vectors
                vec = getproperty(snap, sym)
                if length(size(vec)) == 1
                    f["PartType1/"*name][i, :] = vec
                else
                    f["PartType1/"*name][i, :, :] = vec
                end
            end

            f["Time"][i] = snap.header["Time"]

        end

    end # open f
end

function main()
    din, fout = parse_args()

    check_dir(fout)
    combine_outputs(din, fout)

end # function

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end




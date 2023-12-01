using LilGuys

using Glob

function parse_args()
    if length(ARGS) != 2
        print("expected 2 arguments")
    end
    din, dout = ARGS
    return din, dout
end

function check_dir(dout)
    if ispath(dout)
        println("overwrite $(dout)? (y/N)")
        a = readline()
        if a != "y"
            error("directory exists")
        end
    else
        mkdir(dout)
    end
end

function get_filenames(din)
    path = joinpath(din, "snapshot*.hdf5")
    return glob(path)
end

function main()
    din, dout = parse_args()

    check_dir(dout)
    filenames = get_filenames(din)
    index0 = sort(Snapshot(filenames[1]).index)

    for path in filenames
        _, name = splitdir(path)
        fout = joinpath(dout, name)
        snap = Snapshot(path)
        perm = sortperm(snap.index)
        if index0 != snap.index[perm]
            error("mismatched indicies")
        end
        snap1 = snap[perm]
        LilGuys.write!(fout, snap1)
        println("saving $(fout)")
    end

    paramfile = "parameters-usedvalues"
    in_param = joinpath(din, paramfile)
    cp_param = joinpath(dout, paramfile)
    cp(in_param, cp_param, force=true)
end


main()

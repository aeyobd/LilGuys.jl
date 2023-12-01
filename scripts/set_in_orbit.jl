using ArgParse
using LilGuys

function set_in_orbit(snap, p, v, max_radius=nothing)
    centered = center(snap)
    if max_radius !== nothing
        r = get_r(centered)
        centered = centered[r .< max_radius]
    end
    dp = p ./ R_0 
    dv = v ./ V_0

    centered = shift_snapshot(centered, dp, dv)  # Assuming shift_snapshot is a function you have
    pf, vf = get_most_bound(centered)  # Assuming get_most_bound is a function you have

    println("center at ", pf * R_0)
    println("moving at ", vf * V_0)

    return centered
end

function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "input"
            arg_type=String
            help="input hdf5 file"
        "output"
            arg_type=String
            help="output hdf5 file"
        "-p", "--position"
            arg_type=Float64
            nargs='+'
            help="initial position in kpc"
        "-v", "--velocity"
            arg_type=Float64
            nargs='+'
            help="initial velocity in km/s"
        "--max_radius"
            arg_type=Float64
            default=nothing
            help="clip radius"
    end

    return parse_args(s)
end

function main()
    args = parse_arguments()
    snap = Snapshot(args["input"]) 
    new_snap = set_in_orbit(snap, args["position"], args["velocity"], args["max_radius"])
    save(new_snap, args["output"]) 
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


#!/usr/bin/env julia
#
using ArgParse

using LilGuys 
using HDF5


function get_args()
    s = ArgParseSettings(
        description="Calculate profiles",
    )

    @add_arg_table s begin
        "input"
            help="Input file"
            default="."
        "-o", "--output"
            help="Output file"
            default="profiles.hdf5"
        "-v", "--verbose"
            help="verbose"
            action="store_true"
        "-k", "--skip"
            help="Skip"
            default=10
            arg_type=Int
        "-z", "--zero-centre"
            help="do not use centres"
            action="store_true"

    end

    args = parse_args(s)

    return args
end


function main()
    args = get_args()

    out = Output(args["input"])
    if args["zero-centre"]
        out.x_cen .= zeros(size(out.x_cen))
        out.v_cen .= zeros(size(out.v_cen))
    end

    profiles = LilGuys.ObsProfile3D[]

    snap_idx = eachindex(out)[1:args["skip"]:end]
    for i in snap_idx
        println("computing profile for snapshot $i")
        prof = LilGuys.calc_profile(out[i])
        push!(profiles, prof)
    end


    profs = LilGuys.Profiles3D(snap_idx, profiles)

    LilGuys.save(args["output"], profs)
end


function collect_vector(profiles, field)
    hcat((getfield(p, field) for p in profiles)...)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

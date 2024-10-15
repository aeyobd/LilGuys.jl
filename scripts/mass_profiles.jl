#!/usr/bin/env julia
#
using ArgParse

using LilGuys 
using HDF5

include("bin_args.jl")


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

    s = add_bin_args(s)
    args = parse_args(s)

    return args
end


function main()
    args = get_args()
    bins = bins_from_args(args)

    out = Output(args["input"])
    if args["zero-centre"]
        out.x_cen .= zeros(size(out.x_cen))
        out.v_cen .= zeros(size(out.v_cen))
    end

    profiles = Pair{String, LilGuys.MassProfile3D}[]

    snap_idx = eachindex(out)[1:args["skip"]:end]
    for i in snap_idx
        @info "computing profile for snapshot $i"
        prof = LilGuys.MassProfile3D(out[i], bins=bins)
        push!(profiles, string(i) => prof)
    end


    LilGuys.write_structs_to_hdf5(args["output"], profiles)
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

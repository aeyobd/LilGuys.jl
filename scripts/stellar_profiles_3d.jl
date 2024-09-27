#!/usr/bin/env julia
#
using ArgParse

using LilGuys 
using HDF5


function get_args()
    s = ArgParseSettings(
        description="""
Given the simulation combined output file and a stars file, 
computes the 3D stellar profiles for each snapshot.
""",
    )

    @add_arg_table s begin
        "simulation_output"
            help="Input file"
            required=true
        "starsfile"
            help="Stars file"
            required=true
        "-o", "--output"
            help="Output file"
            default=nothing
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

        # histogram kwargs
        "--bin-method"
            help="Method to calculate the bins"
            default="equal_number"
            arg_type=String
        "--n-bins"
            help="number of bins to use"
            default=nothing
            arg_type=Union{Int, Nothing}

    end

    args = parse_args(s)

    if args["output"] == nothing
        args["output"] = splitext(args["starsfile"])[1] * "_profiles.hdf5"
    end

    return args
end


function main()
    args = get_args()

    weights = LilGuys.read_hdf5_table(args["starsfile"]).probability
    out = Output(args["simulation_output"], weights=weights)

    if args["zero-centre"]
        out.x_cen .= zeros(size(out.x_cen))
        out.v_cen .= zeros(size(out.v_cen))
    end

    profiles = Pair{String, LilGuys.StellarProfile3D}[]

    snap_idx = eachindex(out)[1:args["skip"]:end]

    for i in snap_idx
        @info "computing profile for snapshot $i"
        prof = LilGuys.StellarProfile3D(out[i])
        push!(profiles, string(i) => prof)
    end


    @info "writing to $(args["output"])"
    LilGuys.write_structs_to_hdf5(args["output"], profiles)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

#!/usr/bin/env julia
using ArgParse

using Polyhedra
using LilGuys
import DensityEstimators


include("bin_args.jl")

function get_args()
    s = ArgParseSettings(
        description="Calculate the 2D density profile of a sample of stars in a fits file",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
        help="Input (simulation output) file"
            required=true
        "stars"
            help="Name of the stars table"
            required=true
        "-o", "--output"
            help="Output TOML file"
        "-k", "--skip"
            help="Skip length for snapshots"
            default=10
            arg_type=Int
    end

    s = add_bin_args(s)
    args = parse_args(s)


    if args["output"] === nothing
        samplename = args["input"]
        name = splitext(basename(samplename))[1]
        if dirname(samplename) == ""
            outname = "$(name)_stellar_profiles.toml"
        else
            outname = dirname(samplename) * "/$(name)_stellar_profiles.toml"
        end
        args["output"] = outname
    end

    return args
end



function main()
    args = get_args()
    bins = bins_from_args(args)

    @info "loading sample"
    stars = LilGuys.read_hdf5_table(args["stars"])
    out = Output(args["input"], weights=stars.probability)

    profiles = Pair{String, LilGuys.StellarProfile}[]

    idx = collect(1:args["skip"]:length(out))
    if idx[end] != length(out)
        push!(idx, length(out))
    end

    for i in idx
        @info "processing snapshot $i"
        snap = out[i]

        prof = LilGuys.StellarProfile(snap, bins=bins,
            r_units = "kpc"
          )

        push!(profiles, ("$i" => prof))
    end

    @info "writing data"
    LilGuys.write_structs_to_hdf5(args["output"], profiles)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

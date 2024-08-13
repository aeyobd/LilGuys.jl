#!/usr/bin/env julia

using ArgParse
using LilGuys


function get_args()
    s = ArgParseSettings(
        description="projects a snapshot onto the sky to mock gaia observations",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
            help="Input file"
            required=true
            
        "stars"
            help="Stars file. Expected to be an hdf5 file with columns index, and probabilities with same length as snapshot."
            required=true
        "output"
            help="Output file"

        "-i", "--index"
            help="Index of snapshot to project"
            arg_type=Int
            required=true

        "-m", "--mode"
        help="mode by which to project on sky. May be position (use present-day position to project onto sky), orthoganal (chose an axis and return in physical coordinates) or "
            default="position"
        "-f", "--frame"
        help="frame in which to project on sky. May be ICRS, GSR"
            default="GSR"
    end

    args = parse_args(s)

    return args
end


function get_kwargs(args)
    kwargs = Dict{Symbol, Any}()

    frame = Symbol(args["frame"])
    frame = getproperty(LilGuys, frame)

    kwargs[:SkyFrame] = frame
    return kwargs
end


function main()
    args = get_args()
    kwargs = get_kwargs(args)

    @info "Reading stars"
    stars = LilGuys.load_hdf5_table(args["stars"])

    @assert issorted(stars.index) "stars.index must be sorted"
    @info "Reading snapshot"
    out = Output(args["input"], weights=stars.probabilities)
    snap = out[args["index"]]

    @info "Projecting snapshot onto sky"
    df = LilGuys.to_gaia(snap; kwargs...)

    @info "Writing output"
    LilGuys.write_fits(args["output"], df)

    @info "Done"
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

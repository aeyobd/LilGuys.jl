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
            default="ICRS"
        "-d", "--distance"
            help="distance at which to set snapshot from sun before projecting onto sky"
        "--scale"
            help="file to halo rescaling used.toml"
    end

    args = parse_args(s)

    return args
end


function get_kwargs(args)
    kwargs = Dict{Symbol, Any}()

    frame = Symbol(args["frame"])
    frame = getproperty(LilGuys, frame)

    kwargs[:SkyFrame] = frame

    if args["distance"] !== nothing
        kwargs[:set_to_distance] = parse(Float64, args["distance"])
    end
    return kwargs
end


function main()
    args = get_args()
    kwargs = get_kwargs(args)

    @info "Reading stars"
    stars = LilGuys.read_hdf5_table(args["stars"])

    @assert issorted(stars.index) "stars.index must be sorted"
    @info "Reading snapshot"
    out = Output(args["input"], weights=stars.probability)
    snap = out[args["index"]]

    if args["scale"] != nothing
        scales = TOML.parsefile(args["scale"])
        M_scale = scales["M_scale"]
        v_scale = scales["v_scale"]
        r_scale = scales["r_scale"]

        @assert v_scale^2 == M_scale / r_scale "v_scale^2 must equal M_scale / r_scale"
        snap = LilGuys.rescale(snap, M_scale, r_scale)
    end

    println("snap xcen", snap.x_cen)
    @info "Projecting snapshot onto sky"
    df = LilGuys.to_gaia(snap; kwargs...)
    df_gsr = LilGuys.to_gaia(snap; SkyFrame=LilGuys.GSR)

    df[!, :pmra_gsr] = df_gsr.pmra
    df[!, :pmdec_gsr] = df_gsr.pmdec
    df[!, :radial_velocity_gsr] = df_gsr.radial_velocity

    @info "Writing output"
    LilGuys.write_fits(args["output"], df)

    @info "Done"
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

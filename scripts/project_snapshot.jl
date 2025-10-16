#!/usr/bin/env julia

using ArgParse
using LilGuys
using PyFITS
import TOML

include("script_utils.jl")


function get_args()
    s = ArgParseSettings(
        description="projects a snapshot onto the sky to mock gaia observations",
        version="0.1.1"
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

        "-m", "--mode"
        help="mode by which to project on sky. May be position (use present-day position to project onto sky), orthoganal (chose an axis and return in physical coordinates) or "
            default="position"
        "-f", "--frame"
        help="frame in which to project on sky. May be ICRS, GSR"
            default="ICRS"
        "-d", "--distance"
            help="distance at which to set snapshot from sun before projecting onto sky"
            arg_type=Float64
        "-c", "--centre"
            help="push the centre to the start of the dataframe?"
            action=:store_true
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
        kwargs[:set_to_distance] = args["distance"]
    end

    kwargs[:add_centre] = args["centre"]
    return kwargs
end



function main(args)
    kwargs = get_kwargs(args)
    println(kwargs)

    @info "Reading stars"
    stars = LilGuys.read_hdf5_table(args["stars"])

    @assert issorted(stars.index) "stars.index must be sorted"

    @info "Reading snapshot"
    if args["index"] === nothing
        snap = Snapshot(args["input"])
        if isperm(snap.index)
            idx = snap.index
        else
            idx = invperm(sortperm(snap.index))
        end
        snap.weights = stars.probability[idx]
    else
        out = Output(args["input"], weights=stars.probability)
        snap = out[args["index"]]
    end

    if args["scale"] != nothing
        scales = TOML.parsefile(args["scale"])
        M_scale = scales["M_scale"]
        v_scale = scales["v_scale"]
        r_scale = scales["r_scale"]

        relerr = v_scale^2 / (M_scale / r_scale) - 1
        @assert abs(relerr) < 1e-6 "V scale inconsistent with M & r scale"
        @info "Rescaling snapshot"

        snap = LilGuys.rescale(snap, M_scale, r_scale)
    end

    @info "snap xcen = $(snap.x_cen)"
    @info "Projecting snapshot onto sky"
    kwargs[:filt_wrong_hemisphere] = true
    df = LilGuys.to_gaia(snap; kwargs...)
    kwargs[:SkyFrame] = LilGuys.GSR
    df_gsr = LilGuys.to_gaia(snap; kwargs...)

    df[!, :pmra_gsr] = df_gsr.pmra
    df[!, :pmdec_gsr] = df_gsr.pmdec
    df[!, :radial_velocity_gsr] = df_gsr.radial_velocity

    @info "Writing output"
    write_fits(args["output"], df)

    @info "Done"
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end

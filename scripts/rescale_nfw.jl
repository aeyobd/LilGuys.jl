#!/usr/bin/env julia

using LilGuys
using ArgParse
import TOML

include("script_utils.jl")


function get_args()
    s = ArgParseSettings(description="Rescales a snapshot in mass and radius accordint to the halo parameters.")
    @add_arg_table s begin
        "input"
            help = "Snapshot to rescale."
            required = true
        "--output", "-o"
            help = "file to write snapshot to. Defaults to params"
            default = nothing
            required = true
        "--params-in", "-n"
            help = "TOML file containing the input halo parameters"
            default = nothing
            required = true
        "--params", "-p"
            help = "TOML file containing the desired rescaled halo parameters."
            required = true

        "--max-radius", "-r"
            help = "Truncate particles beyond this radius in the rescaled snapshot."
            default = nothing
            arg_type = Float64
    end

    args = parse_args(s)

    return args
end


function main(args)
    if isfile(args["output"])
        @info "Removing old output file $(args["output"])"
        rm(args["output"], force=true)
    end

    outparams = splitext(args["params"])[1] * "-used.toml"
    rm(outparams, force=true)

    snap = Snapshot(args["input"])

    prof_in = LilGuys.load_profile(args["params-in"])
    params = TOML.parsefile(args["params"])
    prof_out = LilGuys.load_profile(params)

    @assert typeof(prof_in) == typeof(prof_out) "Input and output profiles must be of the same type, got $(typeof(prof_in)) and $(typeof(prof_out)) respectively."

    if (prof_out isa LilGuys.ExpCusp) || (prof_out isa LilGuys.Exp2D)
        m_scale = prof_out.M / prof_in.M
    else
        m_scale = prof_out.M_s / prof_in.M_s
    end

    if prof_out isa LilGuys.Exp2D
        r_scale = prof_out.R_s / prof_in.R_s
    else
        r_scale = prof_out.r_s / prof_in.r_s
    end

    scaled = LilGuys.rescale(snap, m_scale, r_scale)

    if args["max-radius"] !== nothing
        filt = get_r(scaled.positions) .< args["max-radius"]
        @info "Removing $(sum(.!filt)) particles beyond $(args["max-radius"]) kpc"
        scaled = scaled[filt]
    end

    LilGuys.write(args["output"], scaled)

    @info "Writing used scales to $outparams"

    open(outparams, "w") do f
        df = Dict(
            "profile" => params,
            "r_scale" => r_scale,
            "M_scale" => m_scale,
            "v_scale" => sqrt(LilGuys.G*m_scale/r_scale)
           )
        TOML.print(f, df)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end

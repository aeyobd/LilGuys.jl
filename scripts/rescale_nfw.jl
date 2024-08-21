#!/usr/bin/env julia

using LilGuys
using ArgParse
import TOML



function get_args()
    s = ArgParseSettings(description="rescales a snapshot to a new NFW profile.")
    @add_arg_table s begin
        "input"
            help = "snapshot to rescale"
            required = true
        "--output", "-o"
            help = "output snapshot"
            default = nothing
        "--params-in", "-n"
            help = "reads a nfw profile from the file for the input. Defaults to M_s=1, r_s=1"
            default = nothing
        "--params", "-p"
            help = "reads a nfw profile from the file for the output"
            required = true

        "--max-radius", "-r"
            help = "maximum radius to include"
            default = nothing
            arg_type = Float64
    end

    args = parse_args(s)

    if args["output"] === nothing
        args["output"] = splitext(args["params-out"])[1] * ".hdf5"
    end

    return args
end


function main()
    args = get_args()

    snap = Snapshot(args["input"])

    prof_in = LilGuys.load_profile(args["params-in"])
    params = TOML.parsefile(args["params"])
    prof_out = LilGuys.load_profile(params)

    m_scale = prof_out.M_s / prof_in.M_s
    r_scale = prof_out.r_s / prof_in.r_s

    scaled = LilGuys.rescale(snap, m_scale, r_scale)

    if args["max-radius"] !== nothing
        scaled = scaled[get_r(scaled.positions) .< args["max-radius"]]
    end

    LilGuys.save(args["output"], scaled)

    outparams = splitext(args["params"])[1] * "-used.toml"
    println("updating toml at $outparams")
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
    main()
end

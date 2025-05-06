#!/usr/bin/env julia

import Pkg
using Logging, LoggingExtras

using LilGuys
using ArgParse
import TOML


include("script_utils.jl")



function get_args()
    s = ArgParseSettings(description="rescales a snapshot")
    @add_arg_table s begin
        "input"
            help = "snapshot to rescale"
            required = true
        "output"
            help = "output snapshot"
            required = true
        "--mass" , "-m"
            help = "multiplicative mass scale factor"
            arg_type = Float64
        "--radius", "-r"
            help = "multiplicative radius scale factor"
            arg_type = Float64
        "--max-radius"
            help = "truncate particles outside this radius (kpc)"
            default = nothing
            arg_type = Float64
    end

    return parse_args(s)
end


function main(args)
    logfile = splitext(args["output"])[1] * ".log"
    @assert logfile != args["output"]

    logger = TeeLogger(global_logger(), FileLogger(logfile))

    with_logger(logger) do
        Pkg.version()

        if isfile(args["output"])
            rm(args["output"])
        end

        rescale_snapshot(args)
    end
end


function rescale_snapshot(args)
    snap = Snapshot(args["input"])

    r_scale = args["radius"] 
    m_scale = args["mass"]
    scaled = LilGuys.rescale(snap, m_scale, r_scale)

    if args["max-radius"] !== nothing
        scaled = scaled[get_r(scaled.positions) .< args["max-radius"]]
    end

    lguys.save(args["output"], scaled)
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end


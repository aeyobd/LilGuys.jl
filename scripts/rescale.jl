#!/usr/bin/env julia


using LilGuys
using ArgParse
import TOML




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


function main()
    args = get_args()

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
    main()
end


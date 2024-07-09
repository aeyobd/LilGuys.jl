#!/usr/bin/env julia


import LilGuys as lguys
using ArgParse
import TOML


"""
returns a snapshot scaled by the given mass and radius (using code units)
"""
function rescale(snap::lguys.Snapshot, m_scale, r_scale)
    v_scale = sqrt(lguys.G * m_scale / r_scale)

    positions = snap.positions * r_scale
    velocities = snap.velocities * v_scale
    println(typeof(snap.masses))
    masses = snap.masses * m_scale
    println(typeof(masses))
    return lguys.Snapshot(positions, velocities, masses)
end


function get_args()
    s = ArgParseSettings(description="rescales a snapshot")
    @add_arg_table s begin
        "input"
            help = "snapshot to rescale"
            required = true
        "output"
            help = "output snapshot"
            required = true
        "--params", "-p"
            help = "parameter file to read in Ms & Rs"
        "--mass" , "-m"
            help = "scale mass (not mass inside scale radius) in 10^10 Msun"
            arg_type = Float64
        "--radius", "-r"
            help = "scale radius in kpc"
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

    snap = lguys.Snapshot(args["input"])

    if args["params"] !== nothing
        params = TOML.parsefile(args["params"])
        args["mass"] = params["M_s"]
        args["radius"] = params["r_s"]
    end

    r_scale = args["radius"] 
    m_scale = args["mass"] * lguys.A_NFW(1) # zeno scales to M inside r
    scaled = rescale(snap, m_scale, r_scale)

    if args["max-radius"] !== nothing
        scaled = scaled[get_r(scaled.positions) .< args["max-radius"]]
    end

    lguys.save(args["output"], scaled)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


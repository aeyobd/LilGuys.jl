#!/usr/bin/env julia


using LilGuys
using ArgParse
import TOML


function get_args()
    s = ArgParseSettings(description="converts a snapshot to gadget2")
    @add_arg_table s begin
        "input"
            help = "snapshot to rescale"
            required = true
        "output"
            help = "output snapshot"
            required = true
    end

    return parse_args(s)
end


function main()
    args = get_args()

    snap = Snapshot(args["input"])

    if LilGuys.mass_is_fixed(snap.masses)
        m_header = snap.masses[1]
    else
        m_header = 0.0
    end

    snap.header = LilGuys.make_gadget2_header(length(snap), m_header)
    LilGuys.write(args["output"], snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

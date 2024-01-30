#!/usr/bin/env julia

using ArgParse
import LilGuys as lguys

"""
Sets a particle in the orbit given by x_vec_0 and v_vec_0 (kpc and km/s)
"""
function set_in_orbit(snap, x_vec, v_vec, max_radius=nothing)
    centered = lguys.centre(snap)
    if max_radius !== nothing
        r = lguys.calc_r(centered.positions)
        centered = centered[r .< max_radius]
    end

    x0 = x_vec ./ lguys.R0 
    v0 = v_vec ./ lguys.V0

    println("center at ", x_vec)
    println("moving at ", v_vec)
    centered.positions .+= x0
    centered.velocities .+= v0

    return centered
end

function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "input"
            arg_type=String
            help="input hdf5 file"
        "output"
            arg_type=String
            help="output hdf5 file"
        "-p", "--position"
            arg_type=Float64
            nargs='+'
            help="initial position in kpc"
        "-v", "--velocity"
            arg_type=Float64
            nargs='+'
            help="initial velocity in km/s"
        "--max_radius"
            arg_type=Float64
            default=nothing
            help="clip radius"
    end

    return parse_args(s)
end

function main()
    args = parse_arguments()
    snap = lguys.Snapshot(args["input"]) 
    new_snap = set_in_orbit(snap, args["position"], args["velocity"], args["max_radius"])
    lguys.save(args["output"], new_snap) 
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


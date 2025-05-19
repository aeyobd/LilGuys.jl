#!/usr/bin/env julia

using ArgParse
import LilGuys as lguys
using CSV
using DataFrames

include("script_utils.jl")


"""
Sets a particle in the orbit given by x_vec_0 and v_vec_0 (kpc and km/s)
"""
function set_in_orbit(snap, x_vec, v_vec, max_radius=nothing)
    cen = lguys.calc_centre(lguys.SS_State, snap)
    centred = deepcopy(snap)
    centred.positions .-= cen.position
    centred.velocities .-= cen.velocity
    @info "dx cen = $(cen.position)"
    @info "dv cen = $(cen.velocity)"

    if max_radius !== nothing
        r = lguys.calc_r(centred.positions)
        centred = centred[r .< max_radius]
    end

    x0 = x_vec
    v0 = v_vec

    centred.positions .+= x0
    centred.velocities .+= v0

    return centred
end

function get_args()
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
            help="initial position in code units"
        "-v", "--velocity"
            arg_type=Float64
            nargs='+'
            help="initial velocity in code units"
        "-f", "--file"
            arg_type=String
            default=nothing
            help="file to read position and velocity from, overrides -p and -v. Reads the first row as if a csv assuming labeled columns are x, y, .. v_x, v_y, .. etc."
        "--max_radius"
            arg_type=Float64
            default=nothing
            help="clip radius"
    end

    return parse_args(s)
end

function main(args)
    snap = lguys.Snapshot(args["input"]) 
    if args["file"] !== nothing
        df = CSV.read(args["file"], DataFrame)
        row = df[1, :]
        position = [row.x, row.y, row.z]
        velocity = [row.v_x, row.v_y, row.v_z]
        @info "reading from file " * args["file"]
        @info "position = $position"
        @info "velocity = $velocity"
    else
        position = read(file["position"])
        velocity = read(file["velocity"])
    end
    new_snap = set_in_orbit(snap, position, velocity, args["max_radius"])
    lguys.write(args["output"], new_snap) 
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end


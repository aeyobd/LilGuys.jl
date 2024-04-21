#!/usr/bin/env julia

using ArgParse
using JSON

using DataFrames
using CSV
import LilGuys as lguys


function main()
    # Parse command line arguments
    s = ArgParseSettings()

    @add_arg_table s begin
        "input"
            help="Input file"
            default="out/combined.hdf5"
        "-o", "--output"
            help="Output file"
            default="data/centres.csv"
        "-v", "--verbose"
            help="verbose"
            action="store_true"
        "-i", "--maxiter"
            help="Number of iterations"
            arg_type=Int
            default=100
        "-m", "--method"
            help="method to use: shrinking_spheres"
            default="shrinking_spheres"
        "-p", "--percentile"
            help="percentile to keep per round"
            arg_type=Float64
            default=95.
        "-c", "--cut_unbound"
            help="cut unbound particles"
            action="store_true"
        "-f", "--f_min"
            help="minimum fraction of particles"
            default=0.1
            arg_type=Float64
        "-R", "--reinit_state"
            help="do not save previos state"
            action="store_true"
    end


    args = parse_args(s)
    # Read input file
    println("reading snapshots")
    out = lguys.Output(args["input"])

    println("calculating centres")

    kwargs = Dict{Symbol, Any}()
    if args["method"] == "shrinking_spheres"
        statetype = lguys.SS_State
        kwargs[:percen] = args["percentile"]
        kwargs[:f_min] = args["f_min"]
        kwargs[:verbose] = args["verbose"]
    elseif args["method"] == "potential"
        statetype = lguys.StaticState
        kwargs[:method] = "potential"
    elseif args["method"] == "com"
        statetype = lguys.StaticState
        kwargs[:method] = "com"
    end

    kwargs[:reinit_state] = args["reinit_state"]
    cens = lguys.calc_centres(statetype, out; kwargs...)


    x_cen = [cen.position for cen in cens]
    v_cen = [cen.velocity for cen in cens]

    positions = hcat(x_cen...)
    velocities = hcat(v_cen...)

    println("saving centres to ", args["output"])
    df = DataFrame()
    df[!, "t"] = out.times
    df[!, "x"] = positions[1, :]
    df[!, "y"] = positions[2, :]
    df[!, "z"] = positions[3, :]

    df[!, "vx"] = velocities[1, :]
    df[!, "vy"] = velocities[2, :]
    df[!, "vz"] = velocities[3, :]

    CSV.write(args["output"], df)
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

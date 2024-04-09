#!/usr/bin/env julia


using ArgParse
using JSON

import LilGuys as lguys


function main()
    # Parse command line arguments
    s = ArgParseSettings()

    @add_arg_table s begin
        "input"
            help="Input file"
            required=true
        "output"
            help="Output file"
            required=true
        "-v", "--verbose"
            help="verbose"
            action="store_true"
        "-i", "--maxiter"
            help="Number of iterations"
            arg_type=Int
            default=100
        "-p", "--percentile"
            help="percentile to keep per round"
            arg_type=Float64
            default=95.
        "-c", "--cut_unbound"
            help="cut unbound particles"
            action="store_true"
    end


    args = parse_args(s)
    if args["verbose"]
        println(args)
    end
    # Read input file
    input = lguys.Snapshot(args["input"])

    cen = lguys.calc_centre(input, itermax=args["maxiter"], verbose=args["verbose"], percen=args["percentile"])

    new = lguys.copy(input)
    println("shifting by ", cen.x_c, " and ", cen.v_c)
    new.positions .-= cen.x_c
    new.velocities .-= cen.v_c

    if args["cut_unbound"]
        ϵ = lguys.calc_ϵ(new)
        println(ϵ[1:10])
        filt = ϵ .> 0
        println("removing ", sum(ϵ .<= 0), " unbound particles")
        new = new[ϵ .> 0]
        lguys.regenerate_header!(new) # as N has changed
    end

    lguys.save(args["output"], new)
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

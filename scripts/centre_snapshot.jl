#!/usr/bin/env julia


using ArgParse
# using JSON

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
        "-m", "--method"
            help="method to use: shrinking_spheres"
            default="shrinking_spheres"
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

    cen = lguys.calc_centre(statetype, input; kwargs...)

    new = lguys.copy(input)
    println("shifting by ", cen.position, " and ", cen.velocity)
    new.positions .-= cen.position
    new.velocities .-= cen.velocity

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

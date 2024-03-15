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
    end


    args = parse_args(s)
    if args["verbose"]
        println(args)
    end
    # Read input file
    input = lguys.Snapshot(args["input"])

    cen = lguys.ss_centre(input, itermax=args["maxiter"], verbose=args["verbose"], percen=args["percentile"])

    new = lguys.copy(input)
    new.positions .-= cen.x_c
    new.velocities .-= cen.v_c

    Φ = lguys.calc_radial_discrete_Φ(new)
    E = lguys.calc_E_spec_kin(new)
    filt = Φ .+ E .< 0

    new = new[filt]
    lguys.regenerate_header!(new)
    lguys.save(args["output"], new)

end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

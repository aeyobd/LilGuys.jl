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
            default=""
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

    if args["output"] == ""
        println(Dict("x_c" => cen.x_c, "v_c" => cen.v_c))
    else
        open(args["output"], "w") do f
            print(f, Dict("x_c" => cen.x_c, "v_c" => cen.v_c))
        end
    end
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

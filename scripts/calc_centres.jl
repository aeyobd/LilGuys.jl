#!/usr/bin/env julia

using ArgParse

using LilGuys


function main()
    # Parse command line arguments
    s = ArgParseSettings()

    @add_arg_table s begin
        "input"
            help="Input file"
            default="combined.hdf5"
        "-o", "--output"
            help="Output file"
            default="centres.hdf5"
        "-v", "--verbose"
            help="verbose"
            action="store_true"
        "-i", "--maxiter"
            help="Number of iterations"
            arg_type=Int
            default=10
        "-m", "--method"
            help="method to use: ShrinkingSpheres, MostBound, Potential, COM"
            default="ShrinkingSpheres"
        "-q", "--quantile"
            help="fractional percentile to keep per round"
            arg_type=Float64
            default=0.95
        "-c", "--cut_unbound"
            help="cut unbound particles"
            action="store_true"
        "-f", "--f_min"
            help="minimum fraction of particles"
            default=0.001
            arg_type=Float64
        "-R", "--reinit_state"
            help="do not save previos state"
            action="store_true"
        "-k", "--skip"
            help="skip to each nth snapshots"
            arg_type=Int
            default=1

        "--r_max"
            help="maximum radius"
            arg_type=Float64
            default=Inf
        "-a", "--atol"
            help="absolute tolerance"
            arg_type=Float64
            default=1e-3
    end


    args = parse_args(s)
    # Read input file
    println("reading snapshots")
    out = Output(args["input"])

    println("calculating centres")

    kwargs = Dict{Symbol, Any}()
    if args["method"] == "ShrinkingSpheres"
        statetype = LilGuys.SS_State
        kwargs[:r_factor] = args["quantile"] 
        kwargs[:f_min] = args["f_min"]
        kwargs[:verbose] = args["verbose"]
        kwargs[:dx_atol] = args["atol"]
        kwargs[:r_max] = args["r_max"]
    elseif args["method"] == "MostBound"
        statetype = LilGuys.MostBoundState
        kwargs[:percen] = args["percentile"]
        kwargs[:f_min] = args["f_min"]
        kwargs[:verbose] = args["verbose"]
    elseif args["method"] == "Potential"
        statetype = LilGuys.StaticState
        kwargs[:method] = "potential"
    elseif args["method"] == "COM"
        statetype = LilGuys.StaticState
        kwargs[:method] = "com"
    end

    kwargs[:reinit_state] = args["reinit_state"]
    cens = LilGuys.calc_centres(statetype, out; skip=args["skip"], kwargs...)


    x_cen = [cen.position for cen in cens]
    v_cen = [cen.velocity for cen in cens]

    positions = hcat(x_cen...)
    velocities = hcat(v_cen...)

    idx = 1:args["skip"]:length(out)


    df = Dict(
        "snapshots" => idx,
        "times" => out.times,
        "positions" => positions,
        "velocities" => velocities
    )
    write_centres(args["output"], df)
end


function write_centres(filename, df)

end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

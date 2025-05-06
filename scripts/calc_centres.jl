#!/usr/bin/env julia

using ArgParse

using LilGuys

include("script_utils.jl")



function main(args)
    check_overwrite(args)

    @info "reading output"
    out = Output(args["input"])

    @info "calculating centres"
    statetype, kwargs = get_kwargs(args)
    cens = LilGuys.calc_centres(statetype, out; skip=args["skip"], kwargs...)

    println("collecting centres")
    x_cen = [cen.position for cen in cens]
    v_cen = [cen.velocity for cen in cens]

    positions = hcat(x_cen...)
    velocities = hcat(v_cen...)

    idx = collect(1:args["skip"]:length(out))

    df = Dict(
        "snapshots" => idx,
        "times" => out.times[idx],
        "positions" => positions,
        "position_errs" => [cen.position_err for cen in cens],
        "velocities" => velocities,
        "velocity_errs" => [cen.velocity_err for cen in cens],
    )

    @info "writing centres"
    write_centres(args["output"], df)
end


function check_overwrite(args)
    if isfile(args["output"]) && args["interactive"]
        print("File already exists. Overwrite? [y/n] ")
        response = readline()
        if response != "y"
            println("Exiting")
            
            exit()
        end
    end
end


function get_kwargs(args)
    kwargs = Dict{Symbol, Any}()
    if args["method"] == "ShrinkingSpheres"
        statetype = LilGuys.SS_State
        kwargs[:r_factor] = args["quantile"] 
        kwargs[:f_min] = args["f_min"]
        kwargs[:dx_atol] = args["atol"]
        kwargs[:r_max] = args["r_max"]
    elseif args["method"] == "MostBound"
        statetype = LilGuys.MostBoundState
        kwargs[:percen] = args["percentile"]
        kwargs[:f_min] = args["f_min"]
    elseif args["method"] == "Potential"
        statetype = LilGuys.StaticState
        kwargs[:method] = "potential"
    elseif args["method"] == "COM"
        statetype = LilGuys.StaticState
        kwargs[:method] = "com"
    end
    kwargs[:reinit_state] = args["reinit_state"]


    return statetype, kwargs
end


function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "input"
            help="Input file"
            default="combined.hdf5"
        "-o", "--output"
            help="Output file"
            default="centres.hdf5"
        "-I", "--maxiter"
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

        "-i", "--interactive"
            help="prompt before overwriting"
            action="store_true"

        "--r_max"
            help="maximum radius"
            arg_type=Float64
            default=1
        "-a", "--atol"
            help="absolute tolerance"
            arg_type=Float64
            default=1e-3
    end


    args = parse_args(s)
    return args
end


function write_centres(filename, df)
    rm(filename; force=true)

    LilGuys.h5open(filename, "w") do f
        for (key, value) in df
            LilGuys.set_vector!(f, key, value)
        end
    end
end



if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()

    run_script_with_output(main, args)
end

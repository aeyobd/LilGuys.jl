#!/usr/bin/env julia
import Pkg

using Logging, LoggingExtras
using ArgParse

using StatsBase: Histogram, fit, normalize
using HDF5
using CSV, DataFrames
using PyCall

using LilGuys

agama = pyimport("agama")


include("script_utils.jl")


function get_args()
    s = ArgParseSettings(
         description="Animates the dark matter from a (projected) bird's eye view",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
            help="Input (agama .ini) file"

        "-o", "--output"
            help="Output directory for animation frames. Defaults to figures/dm_animation/"
            default="projected_potential.hdf5"
        "-k", "--skip"
            help="Skip length for snapshots"
            default=10
            arg_type=Int
        "-T", "--times"
            help="Times to evaluate projected potential at. May be a list of floats (times to evaluate at), a lguys.Output hdf5 file, or a csv like file with times in the first column."
            nargs='+'
        "-t", "--test"
            help="Test mode: only process the first and last snapshot"
            action="store_true"
        "--limits"
        help="Limits for density calculation. May be a single number (symmetric limits) or four numbers (xmin, xmax, ymin, ymax)"
            nargs='+'
            arg_type=Float64
            default=[200]
        "--x_vec"
            help="x vector for projection"
            nargs=3
            arg_type=Float64
            default=[sind(5), cosd(5), 0]
        "--y_vec"
            help="y vector for projection"
            nargs=3
            arg_type=Float64
            default=[-sind(5), 0, cosd(5)]
        "-n", "--n_bins"
            help="Number of bins (same for both dimensions)"
            default=1001
            arg_type=Int

        "-N", "--n_samples"
            help = "number of samples to use"
            default=1_000_000
            arg_type=Int

        "--r_max"
            help = "maximum radius to sample to"
            default = 1e3
            arg_type = Float64
    end

    args = parse_args(s)
    parse_limits_in_args!(args)
    parse_times_in_args!(args)
    return args
end


"""
Loads in the times to evaluate potential at from the arguments
"""
function parse_times_in_args!(args)
    times_in = args["times"]
    if length(times_in) == 1
        try
            times = [parse(Float64, times_in[1])]
        catch
            filename = times_in[1]
            if endswith(filename, ".hdf5")
                times = Output(filename).times
            else
                # just try it...
                times = CSV.read(filename, DataFrame, header=false)[:, 1] # first column...
            end
        end
    else
        times = [parse(Float64, t) for t in times_in]
    end

    if length(times) == 0 
        error("no times specified")
    end

    args["times"] = times

    return args
end


"""Parses the plot limits from CMD args"""
function parse_limits_in_args!(args)
    if length(args["limits"]) == 1
        args["limits"] = (-args["limits"][1], args["limits"][1], -args["limits"][1], args["limits"][1])
    elseif length(args["limits"]) == 4
        args["limits"] = (args["limits"][1], args["limits"][2], args["limits"][3], args["limits"][4])
    else
        @error "May only specify one or four limits"
    end

    return args["limits"]
end




function project_points(x, y, xybins)
    h1 = fit(Histogram, (x, y), xybins )
    h1 = normalize(h1, mode=:density)
    return h1.weights
end


function sample_potential(pot, time; N, r_max)
    py"""
    import numpy as np
    import agama

    def f(x):
        return np.maximum($pot.density(x, t=$time), 0)
    """

    results = py"agama.sampleNdim(f, $N, np.repeat(-$r_max, 3), np.repeat($r_max, 3))" 

    return results[1]', results[2]

end


function get_xy_samples(pot, time; N, r_max, x_vec, y_vec)

    positions, M = sample_potential(pot, time, N=N, r_max=r_max)

    A = [x_vec y_vec]'
    xy = A * positions

    return xy[1, :], xy[2, :], M
end



function write_single(h5f, h, xbins, ybins, time, x_vec, y_vec, snap)
    attrs = HDF5.attributes(h5f)
    attrs["x_vec"] = x_vec
    attrs["y_vec"] = y_vec
    attrs["time"] = time
    attrs["snapshot"] = snap
    h5f["xbins"] = xbins |> collect
    h5f["ybins"] = ybins |> collect
    h5f["density"] = h
end


function read_and_project(args)
    times = args["times"]

    if isfile(args["output"])
        rm(args["output"])
    end

    pot = agama.Potential(args["input"])

    xbins = LinRange(args["limits"][1], args["limits"][2], args["n_bins"])
    ybins = LinRange(args["limits"][3], args["limits"][4], args["n_bins"])
    bins = (xbins, ybins)
    @assert issorted(xbins) && issorted(ybins)

    if args["test"]
        idx = [1, length(times)]
    else
        idx = 1:args["skip"]:length(times)
    end

    r_max = args["r_max"]

    h5open(args["output"], "w") do f
        for frame in eachindex(idx)
            i = idx[frame]
            @info "processing frame $frame / $(length(idx))"
            x, y, M = get_xy_samples(pot, times[i], r_max=r_max, N=args["n_samples"], x_vec=args["x_vec"], y_vec=args["y_vec"])


            h = project_points(x, y, bins) * M / args["n_samples"]
            dset = "snap$i"

            create_group(f, dset)
            write_single(f[dset], h, xbins, ybins, times[i], args["x_vec"], args["y_vec"], i)
        end
    end
end



function main()
    args = get_args()

    logfile = splitext(args["output"])[1] * ".log"
    @assert logfile != args["output"]

    logger = TeeLogger(global_logger(), FileLogger(logfile))

    with_logger(logger) do
        read_and_project(args)
    end
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

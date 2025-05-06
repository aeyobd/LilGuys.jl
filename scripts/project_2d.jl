#!/usr/bin/env julia
import Pkg
using Logging, LoggingExtras

using ArgParse

using LilGuys

using StatsBase
using LinearAlgebra: cross
using HDF5

include("script_utils.jl")


function get_args()
    s = ArgParseSettings(
         description="Animates the dark matter from a (projected) bird's eye view",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
            help="Input (simulation output) file"
            default="."

        "-o", "--output"
            help="Output directory for animation frames. Defaults to figures/dm_animation/"
            default="projected_densities.hdf5"
        "-s", "--stars"
            help = "path to stars probabilities.hdf5 file"
        "-k", "--skip"
            help="Skip length for snapshots"
            default=10
            arg_type=Int
        "-t", "--test"
            help="Test mode: only process the first and last snapshot"
            action="store_true"
        "--limits"
            help="limits for plot"
            nargs='+'
            arg_type=Float64
            default=[200]
        "--x_vec"
            help="x vector for projection"
            nargs=3
            arg_type=Float64
            default=[0, 1, 0]
        "--y_vec"
            help="y vector for projection"
            nargs=3
            arg_type=Float64
            default=[0, 0, 1]
        "-n", "--n_bins"
            help="number of bins for histogram"
            default=1001
            arg_type=Int
    end

    args = parse_args(s)

    if length(args["limits"]) == 1
        args["limits"] = (-args["limits"][1], args["limits"][1], -args["limits"][1], args["limits"][1])
    elseif length(args["limits"]) == 4
        args["limits"] = (args["limits"][1], args["limits"][2], args["limits"][3], args["limits"][4])
    else
        @error "May only specify one or four limits"
    end

    return args
end



function main(args)
    @info "loading files"

    if args["stars"] != nothing
        df_stars = LilGuys.read_hdf5_table(args["stars"])
        p_idx = df_stars.index
        @assert issorted(p_idx)
        weights = df_stars.probability
    else
        weights = nothing
    end

    out = Output(args["input"], weights=weights)

    xbins = LinRange(args["limits"][1], args["limits"][2], args["n_bins"])
    ybins = LinRange(args["limits"][3], args["limits"][4], args["n_bins"])
    bins = (xbins, ybins)
    @assert issorted(xbins) && issorted(ybins)

    if args["test"]
        idx = [1, length(out)]
    else
        idx = 1:args["skip"]:length(out)
    end

    h5open(args["output"], "w") do f
        for (frame, i) in enumerate(idx)
            print("processing frame $frame / $(length(idx))\r")
            x, y, w = get_xy(out, i; x_vec=args["x_vec"], y_vec=args["y_vec"])
            h = project_points(x, y, w, bins)

            dset = "snap$i"
            create_group(f, dset)
            attrs = HDF5.attributes(f[dset])

            attrs["snapshot"] = i
            attrs["time"] = out.times[i]
            attrs["x_vec"] = args["x_vec"]
            attrs["y_vec"] = args["y_vec"]
            f[dset]["xbins"] = xbins |> collect
            f[dset]["ybins"] = ybins |> collect
            f[dset]["density"] = h

        end
    end
end


function get_xy(out, idx; x_vec=[0, 1, 0], y_vec = [0, 0, 1])
    # shortcut for hdf5
    # -1 since HDF5 is zero indexed
    idx -= 1
    pos = out.h5file["snap$idx/PartType1/Coordinates"][:, :]
    A = [x_vec y_vec]'
    xy = A * pos

    idx = out.h5file["snap$idx/PartType1/ParticleIDs"][:]

    if out.weights != nothing
        if isperm(idx)
            widx = idx
        else
            widx = invperm(sortperm(idx))
        end
        w = out.weights[widx]
    else
        w = nothing
    end

    return xy[1, :], xy[2, :], w
end


function project_points(x, y, w, xybins)
    if w === nothing
        h1 = fit(Histogram, (x, y), xybins )
    else
        h1 = fit(Histogram, (x, y), StatsBase.weights(w), xybins)
    end

    return h1.weights
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end

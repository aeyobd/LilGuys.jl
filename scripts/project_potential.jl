#!/usr/bin/env julia
import Pkg
using Logging, LoggingExtras

using ArgParse

using LilGuys
using StatsBase: midpoints

using LinearAlgebra: cross
using HDF5
using CSV, DataFrames

@info "loading agama"
using PythonCall
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
                times = CSV.read(filename, DataFrame)[:, 1] # first column...
            end
        end
    else
        times = [parse(Float64, t) for t in times_in]
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



function main()
    args = get_args()

    logfile = splitext(args["output"])[1] * ".log"
    @assert logfile != args["output"]

    logger = TeeLogger(global_logger(), FileLogger(logfile))

    with_logger(logger) do
        read_and_project(args)
    end
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

    h5open(args["output"], "w") do f
        if length(idx) == 1
            h = project_agama_potential(pot, bins, time=times[1], x_vec=args["x_vec"], y_vec=args["y_vec"])
            write_single(f, h, xbins, ybins, times[1], args["x_vec"], args["y_vec"], idx[1])
        else
            for frame in eachindex(idx)
                i = idx[frame]
                @info "processing frame $frame / $(length(idx))"
                h = project_agama_potential(pot, bins, time=times[i], x_vec=args["x_vec"], y_vec=args["y_vec"])

                dset = "snap$i"

                create_group(f, dset)
                write_single(f[dset], h, xbins, ybins, times[i], args["x_vec"], args["y_vec"], i)
            end
        end
    end
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

"""
Projects the potential at a given time onto a 2D plane.

"""
function project_agama_potential(potential, bins; x_vec=[0,1,0], y_vec=[0,0,1], time=0)
    xbins, ybins = bins
    xm, ym = midpoints(xbins), midpoints(ybins)
    Nx, Ny = length(xm), length(ym)
    
    alpha, beta, gamma = vectors_to_euler_angles(x_vec, y_vec)

    Σ_disk = Matrix{Float64}(undef, Nx, Ny)

    pos = [repeat(xm, Ny) repeat(ym, inner=Nx)]

    density = potential.projectedDensity(pos, alpha=alpha, beta=beta, gamma=gamma, t=time)

    for i in 1:Ny
        # zero index for python
        idx_i = (i-1)*Nx
        idx_f = idx_i + Nx - 1
        Σ_disk[:, i] .= pyconvert(Vector{Float64}, density[idx_i:idx_f])
    end


    
    return Σ_disk
end


function vectors_to_euler_angles(xhat::Vector{<:Real}, yhat::Vector{<:Real})
    # Compute the zhat vector
    zhat = cross(xhat, yhat)

    # Construct the rotation matrix
    R = hcat(xhat, yhat, zhat)' # transposed so that R * xhat -> x, etc.

    # Extract Euler angles for ZYX convention
    # atan(y, x) or atan(sin, cos)
    alpha = atan(R[3,1], -R[3,2])
    beta = atan(sqrt(R[3,1]^2 + R[3,2]^2), R[3,3])
    gamma = atan(R[1,3], R[2,3])

    return alpha, beta, gamma
end

function project_agama_potential_old(potential, bins; x_vec=[0,1,0], y_vec=[0,0,1], 
        time=0, cutoff=200, cutoffstrength=3)
    # takes a matrix of positions, returns array of densities for each position
    calc_ρ(x) = pyconvert(Vector{Float64}, potential.density(x, t=time))

    xbins, ybins = bins
    xm, ym = midpoints(xbins), midpoints(ybins)
    Nx, Ny = length(xm), length(ym)

    # TODO: cludge...
    zbins = xbins
    zm = midpoints(zbins)
    Nz = length(zm)

    z_vec = cross(x_vec, y_vec)

    Σ_disk = Matrix{Float64}(undef, Nx, Ny)

    for i in 1:Nx # iterate over each x coordinate
        x = fill(xm[i], Ny*Nz)
        # y and z are vectorized so less python calls
        y = repeat(ym, outer=Nz)
        z = repeat(zm, inner=Ny)
        r = @. sqrt(x^2 + y^2 + z^2)
        pos = @. x' * x_vec + y' * y_vec + z' * z_vec
        ρ = calc_ρ(pos')

        if cutoff !== nothing
            ρ = @. ρ * exp(-(r/cutoff)^cutoffstrength)
        end
        ρ = reshape(ρ, Ny, Nz)
        
        Σ = sum(ρ, dims=2) # sum over Z dimension
        Σ_disk[i, :] = Σ
    end

    return Σ_disk
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

#!/usr/bin/env julia
#
import Pkg
using Logging, LoggingExtras

using ArgParse

using LilGuys
using HDF5

import Interpolations: linear_interpolation
import StatsBase: midpoints
import CSV
using DataFrames
import LinearAlgebra: ⋅




include("script_utils.jl")


function get_args()
    s = ArgParseSettings(
         description="resamples the densities in the input file onto the bins and timesteps of the reference file. Useful as preprocessing for animations.",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
            help="Input density evolution file"
        "reference"
            help="reference density file to resample onto"
        "-o", "--output"
            help="Output directory for animation frames. Defaults to figures/dm_animation/"
        "-c", "--centre-file"
            help="Shift the centre of the potential to project by the centre coordinate of this file"
        "-T", "--time-scale"
            help="the relative scale in times as ratio between reference and input times"
            arg_type=Float64
            default=1
    end

    args = parse_args(s)
    return args
end


function main(args)

    if !isnothing(args["centre-file"])
        orbit = get_lmc_orbit(args["centre-file"])
    else
        orbit = nothing
    end

    h5open(args["input"], "r") do f_in
        h5open(args["reference"], "r") do f_ref
            check_consistency(f_in, f_ref)

            snaps_in, times_in = get_snaps_times(f_in)
            snaps_ref = get_snaps(f_ref)

            h5open(args["output"], "w") do f
                for snap in snaps_ref
                    @info "processing $snap"
                    # retrieve relevant values
                    xbins, ybins, Σ = interpolate_density(f_in, f_ref, snap; 
                        snaps_in=snaps_in, times_in=times_in, orbit=orbit,
                        timescale=args["time-scale"])

                    write_new_density(f, f_ref, snap; xbins=xbins, ybins=ybins,
                                      density=Σ)
                end
            end
        end
    end
end


function check_consistency(f_in, f_ref)
    a1 = attrs(f_in["snap1"])
    a2 = attrs(f_ref["snap1"])

    @assert a1["x_vec"] == a2["x_vec"]
    @assert a1["y_vec"] == a2["y_vec"]
end


function interpolate_density(f_in, f_ref, snap; snaps_in, times_in, timescale, orbit)
    s_ref = f_ref[snap]
    time = attrs(s_ref)["time"]
    xbins_new = s_ref["xbins"][:]
    ybins_new = s_ref["ybins"][:]

    # find nearest input snapshots
    i_1, i_2 = get_closest_values(times_in, time / timescale)


    s_1 = f_in[snaps_in[i_1]]
    s_2 = f_in[snaps_in[i_2]]

    t_1 = attrs(s_1)["time"]
    t_2 = attrs(s_2)["time"]

    # interpolate density between snapshots
    x = (time/timescale - t_1) / (t_2 - t_1)

    @info "interpolating $x betweeen $i_1, $i_2"

    Σ1 = s_1["density"][:, :]
    Σ2 = s_2["density"][:, :]
    Σ_m = Σ2 * x .+ Σ1 * (1-x)

    # resample
    xbins, ybins = s_1["xbins"][:], s_1["ybins"][:]

    # apply orbit shift
    if !isnothing(orbit)
        x_hat = attrs(s_1)["x_vec"]
        y_hat = attrs(s_1)["y_vec"]

        centre = [orbit[i](time / timescale) for i in 1:3]
        x_cen = centre ⋅ x_hat 
        y_cen = centre ⋅ y_hat 

        xbins = xbins .+ x_cen
        ybins = ybins .+ y_cen
    end

    # TODO: shift by centre

    Σ_resampled = interpolate_to_grid(Σ_m, xbins, ybins, xbins_new, ybins_new)

    return xbins_new, ybins_new, Σ_resampled
end


function get_lmc_orbit(lmc_file)
    df_lmc = lmc_traj = CSV.read(lmc_file, DataFrame, delim=" ", header = [:time, :x, :y, :z, :v_x, :v_y, :v_z], ignorerepeated=true)

    pos = hcat(df_lmc.x, df_lmc.y, df_lmc.z)'
    vel = hcat(df_lmc.v_x, df_lmc.v_y, df_lmc.v_z)'

    t = df_lmc.time

    return linear_interpolation(t, pos[1, :]), linear_interpolation(t, pos[2, :]), linear_interpolation(t, pos[3, :])
end


function write_new_density(f, f_ref, snap; xbins, ybins, density)

    # make group and write
    create_group(f, snap)

    # set attributes to old values
    a = HDF5.attributes(f[snap])
    for (k, v) in attrs(f_ref[snap])
        a[k] = v
    end

    f[snap]["xbins"] = xbins
    f[snap]["ybins"] = ybins
    f[snap]["density"] = density
end



function get_snaps(f)
    snaps_in = collect(keys(f))
    times_in = [attrs(f[snap])["time"] for snap in snaps_in]

    return snaps_in[sortperm(times_in)]
end


function get_closest_values(times, t)
    idx = searchsortedfirst(times, t)

    return idx-1, idx
end


function get_snaps_times(f)
    snaps_in = collect(keys(f))
    times_in = [attrs(f[snap])["time"] for snap in snaps_in]

    idx = sortperm(times_in)

    return snaps_in[idx], times_in[idx]
end


function interpolate_to_grid(Σ, x, y, xnew, ynew)
    Σ_interp = linear_interpolation((midpoints(x), midpoints(y)), Σ, extrapolation_bc=0.)

    xnew_m = midpoints(xnew)
    ynew_m = midpoints(ynew)

    Σ_new = Matrix{eltype(Σ)}(undef, length(xnew_m), length(ynew_m))


    for j in eachindex(ynew_m)
        Σ_new[:, j] .= Σ_interp.(xnew_m, ynew_m[j])
    end

    return Σ_new
end



if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()

    run_script_with_output(main, args)
end

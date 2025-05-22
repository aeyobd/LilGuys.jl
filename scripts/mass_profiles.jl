#!/usr/bin/env julia
using ArgParse

using LilGuys 
using HDF5
using PyFITS

include("script_utils.jl")


function get_args()
    s = ArgParseSettings(
        description="Calculate profiles",
    )

    @add_arg_table s begin
        "input"
            help="Input file"
            default="."
        "-o", "--output"
            help="Output file"
            default="profiles.hdf5"
        "-v", "--verbose"
            help="verbose"
            action="store_true"
        "-k", "--skip"
            help="Skip"
            default=10
            arg_type=Int
        "-z", "--zero-centre"
            help="do not use centres"
            action="store_true"

    end

    s = add_bin_args(s)
    args = parse_args(s)

    return args
end

function calc_profiles(args)
    outfile_dens = splitext(args["output"])[1] * "_densities.hdf5"
    outfile_scalars = splitext(args["output"])[1] * "_scalars.fits"
    rm(args["output"], force=true)
    rm(outfile_dens, force=true)
    rm(outfile_scalars, force=true)

    bins = bins_from_args(args)

    out = Output(args["input"])

    if args["zero-centre"]
        out.x_cen .= zeros(size(out.x_cen))
        out.v_cen .= zeros(size(out.v_cen))
    end

    snap_idx = eachindex(out)[1:args["skip"]:end]
    profiles = Vector{Pair{String, Union{Missing, LilGuys.MassProfile}}}(undef, length(snap_idx))
    density_profiles = Vector{Pair{String, Union{Missing, LilGuys.DensityProfile}}}(undef, length(snap_idx))
    scalars = Vector{Union{Missing, LilGuys.MassScalars}}(undef, length(snap_idx))

    Threads.@threads for i in eachindex(snap_idx)
        j = snap_idx[i]
        @info "computing profile for snapshot $i"
        try
            prof = LilGuys.MassProfile(out[j], bins=bins)
            dens_prof = LilGuys.DensityProfile(out[j])
            profiles[i] = string(j) => prof
            density_profiles[i] = string(j) => dens_prof
            scalars[i] = MassScalars(out[j], prof)
        catch e
            @warn "failed to compute profile for snapshot $i"
            profiles[i] = string(j) => missing
            density_profiles[i] = string(j) => missing
            scalars[i] = missing
        end
    end

    profiles = filter(x ->!ismissing(last(x)), profiles)
    density_profiles = filter(x ->!ismissing(last(x)), density_profiles)
    scalars = filter(x ->!ismissing(x), scalars)

    LilGuys.write_structs_to_hdf5(args["output"], profiles)
    LilGuys.write_structs_to_hdf5(outfile_dens, density_profiles)
    df_scalars = LilGuys.to_frame(scalars)
    write_fits(outfile_scalars, df_scalars)
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()

    run_script_with_output(calc_profiles, args)
end

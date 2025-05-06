#!/usr/bin/env julia
using ArgParse

using LilGuys 
using HDF5

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
    bins = bins_from_args(args)
    out = Output(args["input"])

    if args["zero-centre"]
        out.x_cen .= zeros(size(out.x_cen))
        out.v_cen .= zeros(size(out.v_cen))
    end

    snap_idx = eachindex(out)[1:args["skip"]:end]
    profiles = Vector{Pair{String, LilGuys.MassProfile3D}}(undef, length(snap_idx))

    Threads.@threads for i in snap_idx
        @info "computing profile for snapshot $i"
        try
            prof = LilGuys.MassProfile3D(out[i], bins=bins)
            profiles[i] = string(i) => prof
        catch e
            @warn "failed to compute profile for snapshot $i"
            @warn e
            profiles[i] = string(i) => nothing
        end
    end

    return profiles
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()

    run_script_with_output(calc_profiles, args)
end

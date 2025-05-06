#!/usr/bin/env julia
#
using ArgParse

using LilGuys 
using HDF5
import TOML

include("script_utils.jl")

SCRIPT_VERSION = "v0.1.0"

function get_args()
    s = ArgParseSettings(
        description="""
Given the simulation combined output file and a stars file, 
computes the 3D stellar profiles for each snapshot.
""",
    )

    @add_arg_table s begin
        "simulation_output"
            help="Input file"
            required=true
        "starsfile"
            help="Stars file"
            required=true
        "-o", "--output"
            help="Output file"
            default=nothing
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
        "--scale"
            help="if set, the file which contains entries for `M_scale` and `r_scale` by which the halo was scaled."

    end

    s = add_bin_args(s)
    args = parse_args(s)

    if args["output"] == nothing
        args["output"] = splitext(args["starsfile"])[1] * "_profiles.hdf5"
    end

    return args
end


function main(args)
    bins = bins_from_args(args)

    weights = LilGuys.read_hdf5_table(args["starsfile"]).probability
    out = Output(args["simulation_output"], weights=weights)
    out_scalar = splitext(args["output"])[1] * "_scalars.hdf5"
    rm(out_scalar, force=true)

    props_file = dirname(args["simulation_output"]) * "/orbital_properties.toml"
    if isfile(props_file)
        idx_peris = TOML.parsefile(props_file)["idx_peris"]
        @assert issorted(idx_peris)
        @info "loaded properties file"
    else
        idx_peris = nothing
    end


    if args["zero-centre"]
        out.x_cen .= zeros(size(out.x_cen))
        out.v_cen .= zeros(size(out.v_cen))
    end

    if args["scale"] != nothing
        scales = TOML.parsefile(args["scale"])
        M_scale = 1 # do not rescale stellar mass...scales["M_scale"]
        M_halo_scale = scales["M_scale"]
        r_scale = scales["r_scale"]
        @info "scaling by M=$M_scale, Mhalo=$M_halo_scale, r=$r_scale"
    else
        M_scale = 1
        M_halo_scale = 1
        r_scale = 1
    end

    profiles = Pair{String, LilGuys.StellarDensity3D}[]
    all_scalars = Pair{String, LilGuys.StellarScalars}[]

    snap_idx = collect(eachindex(out)[1:args["skip"]:end])
    if snap_idx[end] != length(out)
        push!(snap_idx, length(out))
    end

    for i in snap_idx
        @info "computing profile for snapshot $i"
        if (idx_peris != nothing) && (idx_peris[1] <= i)
            idx_last_peri = idx_peris[idx_peris .<= i][end]
            t_peri = out.times[idx_last_peri]
            delta_t = out.times[i] - t_peri
        else
            @info "no peri found"
            delta_t = NaN
        end

        prof = LilGuys.StellarDensity3D(out[i], bins=bins)
        scalars = LilGuys.StellarScalars(out[i], delta_t=delta_t, r_max=1/r_scale)

        @info prof.annotations
        @info "σv original = $(scalars.sigma_v)"

        if args["scale"] != nothing
            prof = LilGuys.scale(prof, r_scale, M_scale)
            scalars = LilGuys.scale(scalars, r_scale, M_scale, M_halo_scale)
            @info "v scaled = $(scalars.sigma_v)"
        end
        push!(profiles, string(i) => prof)
        push!(all_scalars, string(i) => scalars)
    end


    @info "writing to $(args["output"])"
    LilGuys.write_structs_to_hdf5(args["output"], profiles)
    LilGuys.write_structs_to_hdf5(out_scalar, all_scalars)
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end

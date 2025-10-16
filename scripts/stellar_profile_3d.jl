#!/usr/bin/env julia
#
using ArgParse

using LilGuys 
using HDF5
using PyFITS
import TOML

include("script_utils.jl")

SCRIPT_VERSION = "v0.1.0"

function get_args()
    s = ArgParseSettings(
        description="""
Given a snapshot file and a stars file, 
compute the 3D stellar profiles for the snapshot.
""",
    )

    @add_arg_table s begin
        "input"
            help="Input file"
            required=true
        "starsfile"
            help="Stars file"
            required=true
        "output"
            help="Output file"
            required=true
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
        "-u", "--keep-unbound"
            help="keep unbound particles"
            action="store_true"
        "--scale"
            help="if set, the file which contains entries for `M_scale` and `r_scale` by which the halo was scaled."

    end

    s = add_bin_args(s)
    args = parse_args(s)

    return args
end


function main(args)
    bins = bins_from_args(args)

    out_scalar = splitext(args["output"])[1] * "_scalar.toml"
    rm(out_scalar, force=true)

    weights = LilGuys.read_hdf5_table(args["starsfile"]).probability
    snap = Snapshot(args["input"])
    LilGuys.add_stars!(snap, weights)

    props_file = dirname(args["input"]) * "/orbital_properties.toml"
    if isfile(props_file)
        t_peris = TOML.parsefile(props_file)["t_peris"]
        @assert issorted(idx_peris)
        @info "loaded properties file"
    else
        t_peris = nothing
    end


    if args["zero-centre"]
        snap.x_cen .= zeros(3)
        snap.v_cen .= zeros(3)
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

    if (t_peris != nothing) && (t_peris[1] <= snap.time)
        t_peri = maximum(t_peris[t_peris .<= snap.time])
        delta_t = snap.time - t_peri
    else
        @info "no peri found"
        delta_t = NaN
    end

    prof = LilGuys.DensityProfile(snap, snap.weights, bins=bins)
    scalars = LilGuys.StellarScalars(snap, delta_t=delta_t, r_max=1/r_scale)
    if args["keep-unbound"]
        filt_bound = false
    else
        filt_bound = :recursive_1D
    end

    mass_prof = LilGuys.MassProfile(snap, snap.weights, filt_bound=filt_bound)

    @info prof.annotations
    @info "σv original = $(scalars.sigma_v)"

    if args["scale"] != nothing
        prof = LilGuys.scale(prof, r_scale, M_scale, M_halo_scale)
        scalars = LilGuys.scale(scalars, r_scale, M_scale, M_halo_scale)
        @info "v scaled = $(scalars.sigma_v)"
    end

    # save data
    @info "writing to $(args["output"])"
    open(args["output"], "w") do f
        TOML.print(f, LilGuys.struct_to_dict(prof))
    end
    out_mass = splitext(args["output"])[1] * "_mass.toml"
    open(out_mass, "w") do f
        TOML.print(f, LilGuys.struct_to_dict(mass_prof))
    end

    open(out_scalar, "w") do f
        TOML.print(f, LilGuys.struct_to_dict(scalars))
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end

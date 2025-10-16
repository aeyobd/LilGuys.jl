#!/usr/bin/env julia
using ArgParse

using LilGuys 
using HDF5
using PyFITS
import TOML

include("script_utils.jl")


function get_args()
    s = ArgParseSettings(
        description="Calculate profiles",
    )

    @add_arg_table s begin
        "input"
            help="Input file"
        "output"
            help="Output file"
        "-v", "--verbose"
            help="verbose"
            action="store_true"
        "-u", "--keep-unbound"
            help="keep unbound particles"
            action="store_true"
        "-z", "--zero-centre"
            help="do not use centres"
            action="store_true"

    end

    s = add_bin_args(s)
    args = parse_args(s)

    return args
end

function calc_profiles(args)
    outfile_dens = splitext(args["output"])[1] * "_densities.toml"
    outfile_scalars = splitext(args["output"])[1] * "_scalars.toml"
    rm(args["output"], force=true)
    rm(outfile_dens, force=true)
    rm(outfile_scalars, force=true)

    bins = bins_from_args(args)
    snap = Snapshot(args["input"])

    if args["zero-centre"]
        snap.x_cen .= zeros(3)
        snap.v_cen .= zeros(3)
    end

    if args["keep-unbound"]
        filt_bound = false
    else
        filt_bound = :recursive_1D
    end

    prof = LilGuys.MassProfile(snap, bins=bins, filt_bound=filt_bound)
    dens_prof = LilGuys.DensityProfile(snap)
    scalars = MassScalars(snap, prof)

    open(args["output"], "w") do f
        TOML.print(f, LilGuys.struct_to_dict(prof))
    end
    open(outfile_dens, "w") do f
        TOML.print(f, LilGuys.struct_to_dict(dens_prof))
    end

    open(outfile_scalars, "w") do f
        TOML.print(f, LilGuys.struct_to_dict(scalars))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()

    run_script_with_output(calc_profiles, args)
end

#!/usr/bin/env julia

import LilGuys as lguys

using DataFrames
using HDF5

import Statistics: quantile
import TOML
import QuadGK: quadgk

import DensityEstimators

using ArgParse


function get_args()
    s = ArgParseSettings(
        description="""Computes the analytic distribution function
given a snapshot and a halo.
""",
    )

    @add_arg_table s begin
        "energies"
            help="energies file"
            required=true
        "profile"
            help="profile file"
            required=true
        "output"
            help="output file name"
            required=true
        "-n", "--number"
            help="number of radial bins"
            default=10_000
            arg_type=Int
    end

    args = parse_args(s)

    return args
end


function main()
    args = get_args()
    energy_df = lguys.load_hdf5_table(args["energies"])

    prof = lguys.load_profile(args["profile"])
    radii = make_radii_bins(energy_df, args)

    M = sum(energy_df.mass)
    M_scale = M / lguys.calc_M(prof, radii[end]) 
    println("M_scale = $M_scale")

    ρ = M_scale * lguys.calc_ρ.(prof, radii)
    ψ = -M_scale * lguys.calc_Φ.(prof, radii)

    f = lguys.DistributionFunction(ρ, ψ, radii)
    f_ϵ = f.(ψ)

    df_dist = DataFrame(radii=radii, rho=ρ, psi=ψ, f=f_ϵ)

    dist_filename = args["output"] 

    lguys.write_hdf5_table(dist_filename, df_dist; overwrite=true)
end

function make_radii_bins(energy_df, params)
    r = sort(energy_df.radii)

    x_min = log10(r[1])
    x_max = log10(r[end])

    dx_min = log10(r[3]) - log10(r[1])
    dx_max = log10(r[end]) - log10(r[end - 2])

    # pad a little for the derivatives :))
    x_min -= dx_min
    x_max += dx_max

    radii = 10 .^ LinRange(x_min, x_max, params["number"])
    @assert x_max > log10(r[end])
    @assert x_min < log10(r[1])
    return radii
end





if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

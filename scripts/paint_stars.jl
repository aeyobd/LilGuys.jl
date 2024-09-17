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
        description="""Paints stars onto a snapshot.
Based on Rapha's codes which apply Eddingon inversion to determine the 
distribution function of a snapshot.
""",
    )

    @add_arg_table s begin
        "distribution_function"
            help="Name of the distribution function to use"
            required=true
        "energy_file"
            help="File containing the snapshot energies file to paint onto"
            required=true
        "profile"
            help="File containing the profile of the system"
            required=true
        "out"
            help="Output file"
            required=false
    end

    args = parse_args(s)

    if args["out"] == nothing
        args["out"] = splitext(args["profile"])[1]
    end

    return args
end


function main()
    args = get_args()

    profile = lguys.load_profile(args["profile"])
    df_snap = lguys.load_hdf5_table(args["energy_file"])
    df_df = lguys.load_hdf5_table(args["distribution_function"])

    radii = df_df.radii
    ψ = df_df.psi
    print_missing(df_snap.radii, radii, profile)

    ρ = lguys.calc_ρ.(profile, radii)
    f_s = lguys.DistributionFunction(ρ, ψ, radii)


    f_s_e = f_s.(df_df.psi)
    f_dm_e = df_df.f
    prob_e = f_s_e ./ f_dm_e
    if any(prob_e .< 0)
        N_neg = sum(prob_e .< 0)
        @warn "$N_neg/$(length(prob_e)) negative probabilities"

        prob_e[prob_e .< 0] .= 0
    end

    prob = lguys.lerp(ψ, prob_e)

    probs = prob.(df_snap.eps)
    probs = normalize_probabilities(probs)

    df_snap[!, :probability] = probs

    # write outputs
    @info "writing outputs"
    sort!(df_snap, :index)
    lguys.write_hdf5_table(args["out"] * "_stars.hdf5", df_snap; overwrite=true)

    df_df[!, :probability] = prob_e
    df_df[!, :f_s] = f_s_e
    df_df[!, :rho_s] = ρ
    lguys.write_hdf5_table(args["out"] * "_df.hdf5", df_df; overwrite=true)
end


function normalize_probabilities(ps)
    N_neg = sum(ps .< 0)
	@info "$N_neg negative probabilities"
    N_nan = sum(isnan.(ps))
    @info "$N_nan NaN probabilities"
	ps[ps .< 0] .= 0
	ps[isnan.(ps)] .= 0

    if sum(ps) == 0
        error("sum of probabilities is zero")
    end

	ps ./= sum(ps)

    return ps
end


"""
Given the stellar mass function M_s, returns total missing stars
"""
function print_missing(radii, r_bins, profile)
    r_h = lguys.calc_r_h(profile)
    @info "r_h = $r_h"
    @info " $(sum(radii .< r_h)) stars within (3D) half-light radius"

    M_s_tot = lguys.get_M_tot(profile)
    M_s(r) = lguys.calc_M(profile, r) / M_s_tot

	N_s_out = 1 - M_s(r_bins[end])
	N_s_in = M_s(r_bins[1])
	@info "missing $N_s_out stars outside last bin"
	@info "missing $N_s_in stars inside first bin"


    if length(r_bins) < 10
        @warn "only $(length(r_bins)) bins"
        @info "radii: $(r_bins)"
    end
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

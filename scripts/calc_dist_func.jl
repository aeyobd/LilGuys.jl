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
        description="""Computes the distribution function of a snapshot
through eddington inversion.
""",
    )

    @add_arg_table s begin
        "energies"
            help="energies file"
            required=true
        "output"
            help="output file name"
            required=true
    end

    args = parse_args(s)

    return args
end


function main()
    args = get_args()

    energy_df = lguys.load_hdf5_table(args["energies"])

    filt = energy_df.filter

    # radii bins
    radii = energy_df.r[filt]
    r_bins = make_radius_bins(radii, params)
    r_bin_mids = lguys.midpoints(r_bins)


    bins, hist, _ = lguys.histogram(radii, r_bins)
    nu_dm = lguys.calc_ρ_from_hist(bins, hist)
    nu_dm ./= length(radii)
    ϕ = snap_df.phi[snap_df.filter]
    ψ = lguys.lerp(radii, -ϕ).(r_bin_mids)

    df_nu = DataFrame(r=r_bin_mids, nu_dm=nu_dm, nu_s=nu_s, dr=diff(r_bins)/2)

    f_dm = lguys.DistributionFunction(nu_dm, ψ, r_bin_mids)
    @info "writing outputs"

    lguys.write_hdf5_table(params["output_file"] * "_stars.hdf5", snap_df; overwrite=true)

    lguys.write_fits(params["output_file"] * "_density.fits", df_nu)
    lguys.write_fits(params["output_file"] * "_energy.fits", df_E)
end



function calc_phi(snap::lguys.Snapshot)
	radii = lguys.calc_r(snap)
	Φs = lguys.calc_radial_discrete_Φ(radii, snap.masses)
	return Φs
end


function calc_eps(snap, Φs)
	ke = lguys.calc_K_spec(snap)
	ϵs = @. -Φs - ke
	return ϵs
end


function make_filter(snap_df, params)
    ϵs = snap_df.eps
	filt = ϵs .> 0
	return filt
end


"""
given the radii and parameters for stellar profile, returns the
radial bins used in the calculation. 

"""
function make_radius_bins(radii::AbstractVector, params::Dict)

	if !issorted(radii)
		error("radii must be sorted")
	end
	
	r_min = radii[1]
    if keys(params["profile"]) == Set(["KingProfile"])
        prof = lguys.load_profile(params["profile"])
        r_max = prof.R_t
        @info "using King profile with truncation radius $r_max"
        @info "$(sum(radii .< r_max)) particles are within truncation radius"
    else
        r_max = radii[end]
    end
	
    filt = r_min .<= radii
    filt .&= radii .<= r_max
    radii = radii[filt]
	if params["bin_method"] == "equal_width"
		Nr = params["bin_width"]
        log_r = log10(r_min):params["bin_width"]:log10(r_max)
        r_bins = 10 .^ log_r
	elseif params["bin_method"] == "equal_number"
		Nr = params["num_radial_bins"]
		r_bins = quantile(radii, LinRange(0, 1, Nr+1))
	elseif params["bin_method"] == "both"
		r_bins = 10 .^ DensityEstimators.bins_min_width_equal_number(log10.(radii);
		dx_min=params["dr_min"], N_per_bin_min=params["N_per_bin_min"], )
	else
		error("bin method unknown")
	end


	return r_bins
	
end


function make_energy_bins(ψ, params)
	E_max = ψ[1]
	E_min = ψ[end]
	E = LinRange(E_min, E_max, params["num_energy_bins"] + 1)
end


function sample_fs(f_dm, f_s, ψ, params)
    E = make_energy_bins(ψ, params)
    f_dm_e = f_dm.(E)
    f_s_e = f_s.(E)

    if any(f_dm_e .< 0)
        N_neg = sum(f_dm_e .< 0)
        @warn "$N_neg/$(length(f_dm_e)) negative f_dm_e, probabilities may be unreliable"
        E_neg_min = E[findfirst(f_dm_e .< 0)]
        E_neg_max = E[findlast(f_dm_e .< 0)]
        @info "negative f_dm_e span $E_neg_min and $E_neg_max"
    end


    probs = f_s_e ./ f_dm_e
    probs ./= sum(probs .* lguys.gradient(E)) # pdf, dN/dE
    probs[probs .< 0] .= 0

    # set extremes to zero
    E = [0; E;]
    f_dm_e = [0; f_dm_e; ]
    f_s_e = [0; f_s_e; ]
    probs = [0; probs; ]

    return DataFrame(
        E=E,
        f_dm_e=f_dm_e,
        f_s_e=f_s_e,
        probs=probs
    )
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

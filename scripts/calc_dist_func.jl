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
        "-n", "--num-radial-bins"
            help="number of radial bins"
            default=100
            arg_type=Int
        "-b", "--bin-method"
            help="binning method"
            default="equal_number"
        "-p", "--positive"
            help="force DF to be positive"
            action="store_true"
    end

    args = parse_args(s)

    return args
end


function main()
    args = get_args()

    energy_df = lguys.read_hdf5_table(args["energies"])
    sort!(energy_df, :radii)

    filt = energy_df.eps .> 0

    # radii bins
    radii = energy_df.radii[filt]
    r_bins = make_radius_bins(radii, args)
    r_bin_mids = lguys.midpoints(r_bins)
    ψ = energy_df.psi[filt]
    m = energy_df.mass[filt]


    bins, hist, _ = lguys.histogram(radii, r_bins, weights=m)
    ρ = lguys.calc_ρ_from_hist(bins, hist)
    ψ = lguys.lerp(radii, ψ).(r_bin_mids)

    f = lguys.DistributionFunction(ρ, ψ, r_bin_mids; force_positive=args["positive"])
    f_e = f.(ψ)

    df_nu = DataFrame(radii=r_bin_mids, dr=diff(r_bins)/2, rho=ρ, f=f_e, psi=ψ, counts=hist)

    @info "writing outputs"

    lguys.write_hdf5_table(args["output"], df_nu, overwrite=true)
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
    r_max = radii[end]
	
	if params["bin-method"] == "equal_width"
        Nr = params["num-radial-bins"]
        log_r = LinRange(log10(r_min), log10(r_max), Nr+1)
        r_bins = 10 .^ log_r
	elseif params["bin-method"] == "equal_number"
		Nr = params["num-radial-bins"]
		r_bins = quantile(radii, LinRange(0, 1, Nr+1))
	elseif params["bin-method"] == "both"
		r_bins = 10 .^ DensityEstimators.bins_min_width_equal_number(log10.(radii);
		dx_min=params["dr_min"], N_per_bin_min=params["N_per_bin_min"], )
	else
		error("bin method unknown")
	end


	return r_bins
	
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

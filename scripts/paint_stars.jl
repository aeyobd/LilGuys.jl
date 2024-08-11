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
        "input"
            help="Input file"
            required=true
        # "-v", "--verbose"
        #     help="verbose"
        #     action="store_true"
    end

    args = parse_args(s)

    params = load_params(args["input"])
    return params
end


function main()
    params = get_args()
    profile = lguys.load_profile(params)
    snap = load_snap(params)
    snap_df = snap_to_df(snap, params)

    # radii bins
    filt = snap_df.filter
    radii = snap_df.radii[filt]
    r_e = make_radius_bins(radii, params)
    r = lguys.midpoints(r_e)
    print_missing(radii, r_e, profile)

    # density
    _, nu_dm = lguys.calc_ρ_hist(radii, r_e)
    nu_dm ./= length(radii)
    nu_s = max.(lguys.calc_ρ.(profile, r), 0)
    df_nu = DataFrame(r=r, nu_dm=nu_dm, nu_s=nu_s, dr=diff(r_e)/2)

    # distribution functions
    ψ = lguys.lerp(radii, -snap_df.phi[filt]).(r)
    f_dm = lguys.calc_fϵ(nu_dm, ψ, r)
	f_s = lguys.calc_fϵ(nu_s, ψ, r)

    df_E = sample_fs(f_dm, f_s, ψ, params)

    # probabilities
    calc_prob = lguys.lerp(df_E.E, df_E.probs)
    prob = calc_prob.(snap_df.eps[filt])
    prob = normalize_probabilities(prob)
    snap_df[filt, :probability] = prob

    # write outputs
    @info "writing outputs"
    sort!(snap_df, :index)
    lguys.write_hdf5_table(params["output_file"] * "_stars.hdf5", snap_df; overwrite=true)

    lguys.write_fits(params["output_file"] * "_density.fits", df_nu)
    lguys.write_fits(params["output_file"] * "_energy.fits", df_E)
end


"""
loads in the parameterfile and will 
automatically populate the centres file and output files
"""
function load_params(paramname)
    params = TOML.parsefile(paramname); 

    if "output_file" ∉ keys(params)
            params["output_file"] = splitext(paramname)[1]
    end
    
    if "mock_file" ∉ keys(params)
            params["mock_file"] = splitext(paramname)[1] * "_mock_stars.fits"
    end
    
    params
end



function load_snap(params)
    snapname = params["snapshot"]
	@info "loading snap $snapname"
	snap_og = lguys.Snapshot(snapname)

	cen = lguys.calc_centre(lguys.StaticState, snap_og)
	
	snap_og.x_cen = cen.position
	snap_og.v_cen = cen.velocity

    @info "adopting centre $(snap_og.x_cen) $(snap_og.v_cen)"
	
	snap_og
end


function snap_to_df(snap, params)
    snap_df = DataFrame(
        index = snap.index,
        radii = lguys.calc_r(snap),
        phi = calc_phi(snap),
        probability = 0.0,
       )

    snap_df[!, :eps] = calc_eps(snap, snap_df.phi)
    filt = make_filter(snap_df, params)
    snap_df[!, :filter] = filt

    sort!(snap_df, :radii)

    return snap_df
end


function calc_phi(snap)
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
    radii = snap_df.radii

	filt = ϵs .> 0
    if "R_t" in keys(first(params["profile"])[2])
        filt .&= radii .< first(params["profile"])[2]["R_t"]
	end
	return filt
end


"""
given the radii and parameters for stellar profile, returns the
radial bins used in the calculation. 

"""
function make_radius_bins(radii::AbstractVector, params::Dict)
	log_radii = log10.(radii)

	if !issorted(radii)
		error("radii must be sorted")
	end
	
	r_min = radii[1]
	r_max = radii[end]
	
	if params["bin_method"] == "equal_width"
		Nr = params["num_radial_bins"]

		r_e = 10 .^ LinRange(log10.(r_min), log10.(r_max), Nr+1)
	elseif params["bin_method"] == "equal_number"
		Nr = params["num_radial_bins"]
		r_e = quantile(radii, LinRange(0, 1, Nr+1))
	elseif params["bin_method"] == "both"
		r_e = 10 .^ DensityEstimators.bins_min_width_equal_number(log10.(radii);
		dx_min=params["dr_min"], N_per_bin_min=params["N_per_bin_min"], )
	else
		error("bin method unknown")
	end

	return r_e
	
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
    probs = f_s_e ./ f_dm_e
    probs ./= sum(probs .* lguys.gradient(E)) # pdf, dN/dE

    return DataFrame(
        E=E,
        f_dm_e=f_dm_e,
        f_s_e=f_s_e,
        probs=probs
    )
end


function normalize_probabilities(ps)
	@info sum(ps .< 0) " negative probabilities"
	ps[ps .< 0] .= 0
	ps[isnan.(ps)] .= 0
	ps ./= sum(ps)

    return ps
end


"""
Given the stellar mass function M_s, returns total missing stars
"""
function print_missing(radii, r_e, profile)
    # TODO: better determination of these numbers, QuadGK throws a fit if too large a range
    r_min = 1e-5
    r_max = 1000


    r_h = lguys.get_r_h(profile)
    @info " $(sum(radii .< r_h)) stars within (3D) half-light radius"

    ρ_s(r) = lguys.calc_ρ.(profile, r)
	M_s_tot = 4π * quadgk(r-> r^2 * ρ_s(r), r_min, r_max)[1]
	M_s(r) = 4π * quadgk(r-> r^2 * ρ_s(r) / M_s_tot, r_min, r)[1]

	N_s_out = 1 - M_s(r_e[end])
	N_s_in = M_s(r_e[1])
	@info "missing $N_s_out stars outside last bin"
	@info "missing $N_s_in stars inside first bin"
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

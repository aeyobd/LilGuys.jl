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
        "--analytic", "-a"
            help="Pass the name of a profile.toml file to use an analytic density/potential profile instead of the snapshot."
            
    end

    args = parse_args(s)

    params = load_params(args["input"])

    if args["analytic"] != nothing
        params["halo"] = lguys.load_profile(args["analytic"])
    elseif "analytic" ∈ keys(params)
        params["halo"] = lguys.load_profile(params["analytic"])
    end

    return params
end


function main()
    params = get_args()
    profile = lguys.load_profile(params)
    snap = load_snap(params)
    snap_df = snap_to_df(snap, params)

    # radii bins
    radii = snap_df.radii[snap_df.filter]
    r_bins = make_radius_bins(radii, params)
    r_bin_mids = lguys.midpoints(r_bins)
    print_missing(radii, r_bins, profile)
    if length(r_bins) < 10
        @warn "only $(length(r_bins)) bins"
    end

    # density
    if "halo" ∈ keys(params)
        halo = params["halo"]
        @info "using analytic profile"

        nu_dm = lguys.calc_ρ.(halo, r_bin_mids)
        ψ = -lguys.calc_Φ.(halo, r_bin_mids)
        snap_df[:, :phi] = lguys.calc_Φ.(halo, snap_df.radii)
    else
        _, nu_dm = lguys.calc_ρ_hist(radii, r_bins)
        nu_dm ./= length(radii)
        ϕ = snap_df.phi[snap_df.filter]
        ψ = lguys.lerp(radii, -ϕ).(r_bin_mids)
    end

    nu_s = max.(lguys.calc_ρ.(profile, r_bin_mids), 0)
    df_nu = DataFrame(r=r_bin_mids, nu_dm=nu_dm, nu_s=nu_s, dr=diff(r_bins)/2)

    # distribution functions
    f_dm = lguys.DistributionFunction(nu_dm, ψ, r_bin_mids)
	f_s = lguys.DistributionFunction(nu_s, ψ, r_bin_mids)

    df_E = sample_fs(f_dm, f_s, ψ, params)

    # probabilities
    calc_prob = lguys.lerp(df_E.E, df_E.probs)
    prob = calc_prob.(snap_df.eps[snap_df.filter])
    prob[snap_df.eps[snap_df.filter] .< df_E.E[1]] .= 0 # for king profile, stars not bound to inner core can't have probs
    prob = normalize_probabilities(prob)
    snap_df[snap_df.filter, :probability] = prob

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
	log_radii = log10.(radii)

	if !issorted(radii)
		error("radii must be sorted")
	end
	
	r_min = radii[1]
    if keys(params["profile"]) == Set(["KingProfile"])
        r_max = params["profile"]["KingProfile"]["R_t"]
        @info "using King profile with truncation radius $r_max"
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
    N_neg = sum(ps .< 0)
	@info "$N_neg negative probabilities"
	ps[ps .< 0] .= 0
	ps[isnan.(ps)] .= 0
	ps ./= sum(ps)

    return ps
end


"""
Given the stellar mass function M_s, returns total missing stars
"""
function print_missing(radii, r_bins, profile)
    r_h = lguys.calc_r_h(profile)
    @info " $(sum(radii .< r_h)) stars within (3D) half-light radius"

    M_s_tot = lguys.get_M_tot(profile)
    M_s(r) = lguys.calc_M(profile, r) / M_s_tot

	N_s_out = 1 - M_s(r_bins[end])
	N_s_in = M_s(r_bins[1])
	@info "missing $N_s_out stars outside last bin"
	@info "missing $N_s_in stars inside first bin"
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

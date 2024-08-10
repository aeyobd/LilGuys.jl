import LilGuys as lguys

using DataFrames
using HDF5

import StatsBase: percentile
import TOML
import QuadGK: quadgk

using ArgParse


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

    end

    args = parse_args(s)

    return args
end


function main()
    args = get_args()
    params = load_params(paramname)
    profile = load_profile(params)

    return 

    snap_og = load_snap(params)
    snap = filter_and_sort(snap_og)

    # snapshot properties
    ϕ = calc_phi(snap)
    ϵ = calc_eps(snap, ϕ)
    radii = lguys.calc_r(snap)
    M = calc_M_in(radii, snap.masses, r)

    r_e = make_radius_bins(radii, params)
    r = lguys.midpoints(r_e)

    print_missing(radii, r_e, profile)

    _, ν_dm = lguys.calc_ρ_hist(radii, r_e)[2]
    ν_dm ./= length(snap)
    ν_s = max.(lguys.calc_ρ.(profile, r), 0)

    ψ = lguys.lerp(radii, -ϕ).(r)
    f_dm = lguys.calc_fϵ(ν_dm, ψ, r)
	f_s = lguys.calc_fϵ(ν_s, ψ, r)

    df_E = sample_fs(f_dm, f_s, ψ, params)

    calc_prob = lguys.lerp(df_E.E, df_E.probs)
    p = calc_prob.(ϵ)
    p = normalize_probabilities(ps)

    idx_all, ps_all = sort_and_collect(idx_all, idx, ps)
    write_stars(idx_all, ps_all, params)
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
	print(sum(ps .< 0), " negative probabilities")
	ps[ps .< 0] .= 0
	ps[isnan.(ps)] .= 0
	ps ./= sum(ps)

    return ps
end


function sort_and_collect(idx_all, idx, ps)
    idx_excluded = setdiff(idx_all, idx)

    idx_all_sorted = vcat(idx, idx_excluded)
	ps_all = vcat(ps, zeros(length(idx_excluded)))

	_sort = sortperm(idx_all_sorted)
	idx_all_sorted = idx_all_sorted[_sort]
	ps_all = ps_all[_sort]

    @assert idx_all_sorted == sort(snap_og.index) "index does not match snapshot"
    @assert maximum(idx_all_sorted) == length(idx_all) "index does not match length"
    @assert 0 == sum(ps_all[idx_excluded]) "excluded particles have non-zero probability"
    @assert sum(ps_all) == 1" probability sum is not 1"

    return idx_all_sorted, ps_all
end


function load_snap(params)
    snapname = params["snapshot"]
	println("loading snap ", snapname)
	snap_og = lguys.Snapshot(snapname)

	cen = lguys.calc_centre(lguys.StaticState, snap_og)
	
	snap_og.x_cen = cen.position
	snap_og.v_cen = cen.velocity

	println(cen)
	
	snap_og

end


function filter_and_sort(snap_og)
	snap = sort_snapshot(snap_og)

	Φs = calc_phi(snap)
	ϵs = calc_eps(snap, Φs)

	radii = lguys.calc_r(snap)

	filt_snap = make_filter(ϵs, radii, params)
	snap = snap[filt_snap]
end


function sort_snapshot(snap_i)
	radii = lguys.calc_r(snap_i)
	return snap_i[sortperm(radii)]
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


function make_filter(ϵs, radii, params)
	filt = ϵs .> 0
	if "R_t" in keys(params["profile_kwargs"])
		filt .&= radii .< params["profile_kwargs"]["R_t"]
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
		r_e = percentile(radii, LinRange(0, 100, Nr+1))
	elseif params["bin_method"] == "both"
		r_e = 10 .^ Arya.bins_min_width_equal_number(log10.(radii);
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


"""
Given the stellar mass function M_s, returns total missing stars
"""
function print_missing(radii, r_e, profile)
    # TODO: better determination of these numbers, QuadGK throws a fit if too large a range
    r_min = 1e-5
    r_max = 1000


    r_h = lguys.get_r_h(profile)
    println(" $(sum(radii .< r_h)) stars within (3D) half-light radius")

    ρ_s(r) = lguys.calc_ρ.(profile, r)
	M_s_tot = 4π * quadgk(r-> r^2 * ρ_s(r), r_min, r_max)[1]
	M_s(r) = 4π * quadgk(r-> r^2 * ρ_s(r) / M_s_tot, r_min, r)[1]

	N_s_out = 1 - M_s(r_e[end])
	N_s_in = M_s(r_e[1])
	println("missing $N_s_out stars outside")
	println("missing $N_s_in stars inside")
end


"""
Writes the stars to an hdf5 file
"""
function write_stars()
	outname = params["output_file"]
	if isfile(outname)
		if overwrite
			rm(outname)
		else
			println("file already exists")
			return
		end
	end


	h5write(outname, "/index", idx_all)
	h5write(outname, "/probabilities", ps_all)
	println("probabilities written to $outname")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end



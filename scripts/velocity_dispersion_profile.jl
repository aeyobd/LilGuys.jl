#!/usr/bin/env julia
using ArgParse

import Polyhedra
using LilGuys
using PyFITS
import TOML
using OrderedCollections

import StatsBase

include("script_utils.jl")

function get_args()
    s = ArgParseSettings(
        description="Calculate velocity dispersion profile with elliptical radiusof a dataframe",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
            help="Input fits file"
            required=true
        "output"
            help="Output TOML file"
        "-s", "--simulation"
            help="simulation flags"
            action="store_true"
        "--ellipticity"
            help="Ellipticity of the profile"
            default=0
            arg_type=Float64
        "--PA"
            help="Position angle of the profile"
            default=0
            arg_type=Float64
        "--centre-method"
            help="""Method to calculate the centre of the profile. May be
                mean, weighted, weighted3, or first."""
            default="mean"
            arg_type=String
        "--mass-column"
            help="Column in the fits file with the mass of the stars"
            default=nothing
            arg_type=String
    end

    s = add_bin_args(s)
    args = parse_args(s)


    if args["output"] === nothing
        samplename = args["input"]
        name = splitext(basename(samplename))[1]
        if dirname(samplename) == ""
            outname = "$(name)_sigma_v_profile.toml"
        else
            outname = dirname(samplename) * "/$(name)_sigma_v_profile.toml"
        end
        args["output"] = outname
    end

    if args["simulation"]
        args["mass-column"] = "weights"
        args["centre-method"] = "weighted3"
    end


    return args
end



function main()
    args = get_args()
    bins = bins_from_args(args)

    @info "loading sample"
    sample = read_fits(args["input"])

    if args["mass-column"] === nothing
        weights = ones(size(sample, 1))
    else
        @info "setting weights"
        weights = sample[:, args["mass-column"]]
    end

    @info "calculating centre"
    ra0, dec0 = calc_centre(sample, weights, args)
    @info "centre: $(ra0) $(dec0)"

    @info "calculating profile"
    r_ell = LilGuys.calc_R_ell_sky(sample.ra, sample.dec,
        args["ellipticity"], args["PA"], centre=(ra0, dec0))

    log_R_bins = bins(log10.(r_ell), weights)
    log_R = LilGuys.midpoints(log_R_bins)
    R_bins = exp10.(log_R_bins)

    σs, err, Ns = calc_sigma_v_profile(r_ell, sample.radial_velocity, weights; R_bins=R_bins)

    df_out = OrderedDict(
        "R_bins" => R_bins,
        "sigma_v" => σs,
        "sigma_v_err" => err,
        "counts" => Ns,
        "R" => exp10.(log_R),
       )

    open(args["output"], "w") do f
        TOML.print(f, df_out)
    end

    @info "wrote data to $(abspath(args["output"]))"
end




function calc_sigma_v_profile(R_ell, radial_velocity, weights; R_bins)
    filt = @. isfinite(R_ell)

    r = R_ell[filt]
	w = weights[filt]
	v = radial_velocity[filt]

    N_bins = length(R_bins) - 1
	σs = Vector{Float64}(undef, N_bins)
	err = Vector{Float64}(undef, N_bins)
	Ns = Vector{Float64}(undef, N_bins)
	
	for i in 1:N_bins
		filt = R_bins[i] .<= r .< R_bins[i+1]
        σs[i] = LilGuys.std(v[filt], StatsBase.weights(w[filt]))
		N = sum(filt)
        err[i] = NaN
		Ns[i] = N
	end

    return σs, err, Ns
end

function calc_centre(sample, weights, args)
    if args["centre-method"] == "mean"
        ra0, dec0 = LilGuys.calc_centre2D(sample.ra, sample.dec, "mean")
    elseif args["centre-method"] == "weighted"
        ra0, dec0 = LilGuys.calc_centre2D(sample.ra, sample.dec, "mean", weights)
    elseif args["centre-method"] == "weighted3"
        ra0, dec0 = LilGuys.calc_centre2D(sample.ra, sample.dec, "mean", weights .^ 3)
    elseif args["centre-method"] == "first"
        ra0 = sample.ra[1]
        dec0 = sample.dec[1]
    else
        throw(ArgumentError("Invalid centre method $(args["centre-method"])"))
    end

    return ra0, dec0
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

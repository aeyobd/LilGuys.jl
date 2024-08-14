#!/usr/bin/env julia
using ArgParse

using Polyhedra
using LilGuys
import DensityEstimators


function get_args()
    s = ArgParseSettings(
        description="Calculate the 2D density profile of a sample of stars in a fits file",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
            help="Input fits file"
            required=true
        "output"
            help="Output TOML file"
        "--r_centre"
            help="Centre of the profile"
            default=5
        "--N-min"
            help="Minimum number of stars per bin"
            default=20
            arg_type=Int
        "--dlogr-min"
            help="Minimum bin width"
            default=0.05
            arg_type=Float64
        "--ellipticity"
            help="Ellipticity of the profile"
            default=0
            arg_type=Float64
        "--PA"
            help="Position angle of the profile"
            default=NaN
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
        "--normalization"
            help="Normalization of the profile"
            default="central"
            arg_type=String
        "--r-centre"
            help="Radius past which to normalize"
            default=5
            arg_type=Float64
    end

    args = parse_args(s)

    if args["output"] === nothing
        samplename = args["input"]
        name = splitext(basename(samplename))[1]
        if dirname(samplename) == ""
            outname = "$(name)_profile.toml"
        else
            outname = dirname(samplename) * "/$(name)_profile.toml"
        end
        args["output"] = outname
    end


    if args["ellipticity"] > 0 && isnan(args["PA"])
        throw(ArgumentError("Ellipticity is set but PA is not"))
    elseif args["ellipticity"] == 0 && isnan(args["PA"])
        args["PA"] = 0
    end

    return args
end



function main()
    args = get_args()

    @info "loading sample"
    sample = LilGuys.load_fits(args["input"])

    if args["mass-column"] === nothing
        weights = ones(size(sample, 1))
    else
        @info "setting weights"
        weights = sample[:, args["mass-column"]]
        println(maximum(weights))
    end

    @info "calculating centre"
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


    @info "centre: $(ra0) $(dec0)"

    @info "calculating profile"
    r_ell = LilGuys.calc_r_ell_sky(sample.ra, sample.dec,
        args["ellipticity"], args["PA"], centre=(ra0, dec0))

    r_ell_max = LilGuys.calc_r_max(sample.ra, sample.dec, 
        args["ellipticity"], args["PA"], centre=(ra0, dec0))

    filt = r_ell .< r_ell_max

    bins =  DensityEstimators.bins_min_width_equal_number(log10.(r_ell[filt]), 
        N_per_bin_min=args["N-min"], dx_min=args["dlogr-min"])

    profile = LilGuys.calc_properties(r_ell[filt], bins=bins, weights=weights[filt],
        normalization=Symbol(args["normalization"]), r_centre=args["r-centre"])


    @info "writing data"
    open(args["output"], "w") do f
		print(f, profile)
	end

    @info "wrote data to ", abspath(args["output"])
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
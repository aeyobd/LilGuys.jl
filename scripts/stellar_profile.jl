#!/usr/bin/env julia
using ArgParse

using Polyhedra
using LilGuys
using PyFITS

include("script_utils.jl")

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
        "-s", "--simulation"
            help="simulation flags"
            action="store_true"
        "--r_centre"
            help="Centre of the profile"
            default=5
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
            default="none"
            arg_type=String
        "--r-centre"
            help="Radius past which to normalize"
            default=5
            arg_type=Float64
        "--rv-max-radius"
            help="Maximum radius for radial velocity calculation in arcminutes"
            default=60
            arg_type=Float64
    end

    s = add_bin_args(s)
    args = parse_args(s)


    if args["simulation"]
        args["mass-column"] = "weights"
        args["centre-method"] = "weighted3"
    end

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
    r_ell = LilGuys.calc_r_ell_sky(sample.ra, sample.dec,
        args["ellipticity"], args["PA"], centre=(ra0, dec0))

    filt = .!isnan.(r_ell)

    if args["simulation"]
        distance = LilGuys.mean(sample.distance, weights)
        @info "distance: $(distance)"

        rv_filt = sample.r_ell .< args["rv-max-radius"]

        σv = LilGuys.std(sample.radial_velocity[rv_filt], weights[rv_filt])
        @info "σv: $(σv)"
    else
        distance = NaN
        σv = NaN
    end

    profile = LilGuys.StellarProfile(r_ell[filt], 
        bins=bins, 
        weights=weights[filt], 
        normalization=Symbol(args["normalization"]), 
        r_centre=args["r-centre"],
        distance=distance,
        sigma_v=σv,
    )


    @info "writing data"
    open(args["output"], "w") do f
		print(f, profile)
	end

    @info "wrote data to ", abspath(args["output"])
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

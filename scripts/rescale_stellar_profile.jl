#!/usr/bin/env julia
#
using LilGuys
using ArgParse


function get_args()
    s = ArgParseSettings(description="rescales a stellar profile")
    @add_arg_table s begin
        "input"
            help = "stellar profile to rescale"
            required = true
        "output"
            help = "output stellar profile"
            required = true
        "--mass" , "-m"
        help = "mass scale (as factor)"
            arg_type = Float64
        "--radius", "-r"
            help = "radius scale (as factor)"
            arg_type = Float64
        "--radius-units", "-u"
            help = "units of resulting radii"
            arg_type = String
    end

end


function main()
    args = get_args()
    prof = StellarProfile(args["input"])

    r_scale = args["radius"]
    m_scale = args["mass"]
    scaled = LilGuys.rescale_stellar_profile(prof, r_scale, m_scale)

    scaled.r_units = args["radius-units"]

    lguys.save(args["output"], scaled)
end

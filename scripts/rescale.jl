using LilGuys
using ArgParse


function rescale(snap::Snapshot, m_scale, r_scale)
    scaled = copy(snap)
    scaled.pos .*= r_scale
    scaled.vel .*= sqrt(LilGuys.G * m_scale / r_scale)
    scaled.m *= m_scale
    return scaled
end


function get_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--mass" , "-m"
            help = "scale mass in 1e10 Msun"
            arg_type = F
        "--radius", "-r"
            help = "scale radius in kpc"
            arg_type = F
        "--max-radius"
            help = "truncate particles outside this radius (kpc)"
            default = nothing
            arg_type = F
        "input"
            help = "snapshot to rescale"
            required = true
        "output"
            help = "output snapshot"
            required = true
    end

    return parse_args(s)
end


function main()
    args = get_args()

    snap = Snapshot(args["input"])

    r_scale = args["radius"] / LilGuys.R0
    m_scale = args["mass"]/ LilGuys.M0
    v_scale = sqrt(LilGuys.G * m_scale / r_scale)

    scaled = copy(snap)
    scaled.pos .*= r_scale
    scaled.vel .*= v_scale
    scaled.m *= m_scale

    if args["max-radius"] !== nothing
        scaled = scaled[get_r(scaled.pos) .< args["max-radius"]]
    end

    write!(args["output"], scaled)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


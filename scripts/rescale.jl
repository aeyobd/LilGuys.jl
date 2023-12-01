using LilGuys
using ArgParse


function get_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--mass" , "-m"
            help = "scale mass in 1e10 Msun"
        "--radius", "-r"
            help = "scale radius in kpc"
        "--max-radius"
            help = "truncate particles outside this radius (kpc)"
            default = nothing
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

    r_scale = args["radius"] / LilGuys.R_0
    m_scale = args["mass"]/ LilGuys.M_0
    v_scale = sqrt(LilGuys.G * m_scale / LilGuys.R_0)

    scaled = copy(snap)
    scaled.pos .*= r_scale
    scaled.vel .*= v_scale
    scaled.m *= m_scale

    if args.max_radius !== nothing
        scaled = scaled[get_r(scaled.pos) .< args.max_radius]
    end

    save(scaled, args["output"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


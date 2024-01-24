import LilGuys as lguys
using ArgParse


function rescale(snap::lguys.Snapshot, m_scale, r_scale)
    v_scale = sqrt(lguys.G * m_scale / r_scale)

    scaled = lguys.copy(snap)
    scaled.positions .*= r_scale
    scaled.velocities .*= v_scale
    if scaled.masses isa lguys.ConstVector
        scaled.masses *= m_scale
    else
        scaled.masses .*= m_scale
    end
    return scaled
end


function get_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--mass" , "-m"
            help = "scale mass in 1e10 Msun"
            arg_type = Float64
        "--radius", "-r"
            help = "scale radius in kpc"
            arg_type = Float64
        "--max-radius"
            help = "truncate particles outside this radius (kpc)"
            default = nothing
            arg_type = Float64
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

    snap = lguys.Snapshot(args["input"])

    r_scale = args["radius"] / lguys.R0
    m_scale = args["mass"]/ lguys.M0
    scaled = rescale(snap, m_scale, r_scale)

    if args["max-radius"] !== nothing
        scaled = scaled[get_r(scaled.positions) .< args["max-radius"]]
    end

    lguys.save(args["output"], scaled)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


#!/usr/bin/env julia
#
using ArgParse

using LilGuys 
using HDF5

function main()
    args = get_args()

    out = Output(args["input"])

    profiles = LilGuys.ObsProfile3D[]

    snap_idx = eachindex(out)[1:args["skip"]:end]
    for i in snap_idx
        println("computing profile for snapshot $i")
        prof = LilGuys.calc_profile(out[i])
        push!(profiles, prof)
    end

    # scalars
    E = [p.E for p in profiles]
    W = [p.W for p in profiles]
    K = [p.K for p in profiles]
    v_circ_max = [p.v_circ_max for p in profiles]
    r_circ_max = [p.r_circ_max for p in profiles]
    N_bound = [p.N_bound for p in profiles]

    # vectors
    log_r = collect_vector(profiles, :log_r)
    rho = collect_vector(profiles, :ρ)
    v_circ = collect_vector(profiles, :v_circ)
    M_in = collect_vector(profiles, :M_in)

    h5open(args["output"], "w") do f
        f["snapshots"] = collect(snap_idx)
        f["E"] = E
        f["W"] = W
        f["K"] = K
        f["v_circ_max"] = v_circ_max
        f["r_circ_max"] = r_circ_max
        f["N_bound"] = N_bound

        f["log_r"] = log_r

    end

end


function collect_vector(profiles, field)
    hcat((getfield(p, field) for p in profiles)...)
end

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



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

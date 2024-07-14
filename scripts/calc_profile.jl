#!/usr/bin/env julia
#
using ArgParse
import TOML

using LilGuys 

function main()
    args = get_args()
    snap = load_snap(args["input"], args["snap"], args["centres"])
    df = global_properties(snap)
end


function global_properties(snap)

    E_kin = lguys.calc_K_tot(snap)
    E_pot = lguys.calc_W_tot(snap)
    E_tot = E_kin + E_pot


    fit = lguys.fit_v_circ(r, v_circ)

    df = Dict(
        "E_kin" => E_kin,
        "E_pot" => E_pot,
        "E_tot" => E_tot,
        "v_circ_max" => fit[:v_circ_max],
        "r_circ_max" => fit[:r_circ_max],
    )
    return df
end


function profiles(snap)
    r, v_circ = lguys.calc_v_circ(snap)


    df = Dict(
        "r_circ" => r,
        "v_circ" => v_circ
    )
end


function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "input"
            help="Input file"
            default="combined.hdf5"
        "-o", "--output"
            help="Output file"
            default="centres.csv"
        "-v", "--verbose"
            help="verbose"
            action="store_true"
        "-c" "--centres"
            help="Centres file"
            default="centres.hdf5"
        "-s", "--snap"
            help="Snap index"
            default=1
            arg_type=Int

    end

    args = parse_args(s)

    return args
end


function load_snap(filename, snap, centres)
    out = Output(filename)

    LilGuys.h5open(centres, "r") do f
        if LilGuys.get_vector(f, "snapshots")[end] != length(out)
            error("out index does not match")
        end

        out.x_cen .= LilGuys.get_vector(f, "positions")
        out.v_cen .= LilGuys.get_vector(f, "velocities")
    end

    return out[snap]
end




if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

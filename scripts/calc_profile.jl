#!/usr/bin/env julia
#
using ArgParse
import TOML

using LilGuys 

function main()
    args = get_args()
    snap = load_snap(args["input"], args["snap"])
    df = global_properties(snap)
    df = profiles(snap)
end


function global_properties(snap)

    E_kin = LilGuys.calc_K_tot(snap)
    E_pot = LilGuys.calc_W_tot(snap)
    E_tot = E_kin + E_pot


    fit = LilGuys.fit_v_r_circ_max(snap)

    df = Dict(
        "E_kin" => E_kin,
        "E_pot" => E_pot,
        "E_tot" => E_tot,
        "v_circ_max" => fit[:v_circ_max],
        "r_circ_max" => fit[:r_circ_max],
    )

    println(df)
    return df
end


function profiles(snap; N_bins=100)
    radii = LilGuys.calc_radii(snap)
    idx_sorted = sortperm(radii)
    M = cumsum(snap.masses)[idx]

    r_rho, rho = LilGuys.calc_ρ_hist(snap, N_bins)

    df = Dict(
        "r_circ" => r_circ,
        "v_circ" => v_circ,
        "r_rho" => r_rho,
        "rho" => rho,
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
        "-v", "--verbose"
            help="verbose"
            action="store_true"
        "-s", "--snap"
            help="Snap index"
            default=1
            arg_type=Int

    end

    args = parse_args(s)

    return args
end


function load_snap(filename, snap)
    out = Output(filename)
    return out[snap]
end




if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

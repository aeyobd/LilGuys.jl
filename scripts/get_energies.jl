#!/usr/bin/env julia

import LilGuys as lguys

using DataFrames
using HDF5

import Statistics: quantile
import TOML
import QuadGK: quadgk

import DensityEstimators

using ArgParse


function get_args()
    s = ArgParseSettings(
        description="""Calculates the energies of a snapshot, 
outputing a table of the snapshot index, binding energy, potential, and
radius for each particle.
""",
    )

    @add_arg_table s begin
        "snapshot"
            help="snapshot to use"
            required=true
        "output"
            help="output file name"
            required=true
    end

    args = parse_args(s)

    return args
end


function main()
    args = get_args()
    snap = lguys.Snapshot(args["snapshot"])

	cen = lguys.calc_centre(lguys.StaticState, snap)
    # radii bins
    radii = lguys.calc_r(snap, cen.position)
    masses = snap.masses
    Φs = lguys.calc_radial_discrete_Φ(radii, masses)
    vs = lguys.calc_r(snap.velocities, cen.velocity)
    ϵs = @. -Φs - 1/2*vs^2

    snap_df = DataFrame(
        index = snap.index,
        radii = radii,
        mass = masses,
        psi = -Φs,
        psi_nbody = -snap.Φs,
        eps = ϵs,
        filt = ϵs .> 0,
    )

    lguys.write_hdf5_table(args["output"] , snap_df; overwrite=true)

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

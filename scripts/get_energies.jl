#!/usr/bin/env julia

import LilGuys as lguys

using DataFrames
using HDF5

import Statistics: quantile
import TOML
import QuadGK: quadgk

import DensityEstimators

using ArgParse

include("script_utils.jl")


function get_args()
    s = ArgParseSettings(
        description="""Calculates the energies of a snapshot, 
outputing a table of the snapshot index, binding energy, potential, and
radius for each particle.
Note, best practice is to pull the snapshot from a combined.hdf5 file with centres.
""",
    )

    @add_arg_table s begin
        "snapshot"
            help="snapshot to use"
            required=true
        "output"
            help="file to write the energies to"
            required=true
    end

    args = parse_args(s)

    return args
end


function main(args)
    snap = lguys.Snapshot(args["snapshot"])

    cen = lguys.calc_centre(lguys.StaticState, snap)
    # radii bins
    radii = lguys.radii(snap, cen.position)
    masses = snap.masses
    Φs = lguys.potential_spherical_discrete(radii, masses)
    vs = lguys.radii(snap.velocities, cen.velocity)
    ϵs = @. -Φs - 1/2*vs^2
    Φs_nbody = isnothing(snap.potential) ? fill(NaN, length(Φs)) : snap.potential

    snap_df = DataFrame(
        index = snap.index,
        radii = radii,
        mass = masses,
        psi = -Φs,
        psi_nbody = Φs_nbody,
        eps = ϵs,
        filt = ϵs .> 0,
    )

    lguys.write_hdf5_table(args["output"] , snap_df; overwrite=true)
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    run_script_with_output(main, args)
end

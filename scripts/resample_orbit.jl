#!/usr/bin/env julia

using ArgParse
using CSV, DataFrames
using LilGuys

function get_args()
    s = ArgParseSettings(description="resamples a CSV orbit")
    @add_arg_table s begin
        "input"
            help = "Path to the LMC trajectory file"
            required = true
        "--output", "-o"
            help = "Output file"
            required = true
        "--times", "-t"
            help = "path to an outputfile (if .hdf5) a orbit file (if .CSV) or a list of tiles"
            required = true
        "--time-scale"
            arg_type = Float64
            default = 1.0
            help = "Time scale factor"
    end

    return parse_args(s)
end

function main()
    args = get_args()

    # Parse inputs
    filename = args["input"]
    outputname = args["output"]

    # Read LMC trajectory file
    lmc_traj = CSV.read(filename, DataFrame, delim=" ", header=[:time, :x, :y, :z, :v_x, :v_y, :v_z])

    # Constants
    V_T2GYR = 0.97779

    # Scale time
    lmc_traj.time .= lmc_traj.time .* V_T2GYR / T2GYR

    timesin  = args["times"]
    if timesin isa AbstractString
        if endswith(timesin, ".csv")
            times = CSV.read(timesin, DataFrame).t
        elseif endswith(timesin, ".hdf5")
            times = Output(timesin).times
        end
    else
        times = parse.(Float64, timesin)
    end 

    times = times .* args["time-scale"]

    # Interpolate trajectories
    lmc_x = LilGuys.lerp(lmc_traj.time, lmc_traj.x).(times)
    lmc_y = LilGuys.lerp(lmc_traj.time, lmc_traj.y).(times)
    lmc_z = LilGuys.lerp(lmc_traj.time, lmc_traj.z).(times)

    lmc_vx = LilGuys.lerp(lmc_traj.time, lmc_traj.v_x).(times)
    lmc_vy = LilGuys.lerp(lmc_traj.time, lmc_traj.v_y).(times)
    lmc_vz = LilGuys.lerp(lmc_traj.time, lmc_traj.v_z).(times)

    # Create DataFrame
    lmc = DataFrame(
        :time => times,
        :x => lmc_x,
        :y => lmc_y,
        :z => lmc_z,
        :v_x => lmc_vx ./ V2KMS,
        :v_y => lmc_vy ./ V2KMS,
        :v_z => lmc_vz ./ V2KMS,
    )

    # Write output
    CSV.write(outputname, lmc)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

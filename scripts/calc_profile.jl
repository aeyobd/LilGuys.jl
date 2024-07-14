
using ArgParse
import TOML

using LilGuys 
using DataFrames, CSV

function main()
    args = parse_args()
    snap = Snapshot(args["input"])
    cen = 
end


function parse_args()
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
    end

    args = parse_args(s)

    return args
end






if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

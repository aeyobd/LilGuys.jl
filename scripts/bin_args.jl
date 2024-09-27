using LilGuys: bins_equal_width, bins_equal_number, bins_both
using ArgParse


"""
    bins_from_args

given `args` with options as in add_bin_args below, returns 
the bins function to pass to a histogram.

"""
function bins_from_args(args; logarithmic=false)
    binmethod = args["bin-method"]

    local f_bins
    if binmethod == "equal-number"
        f_bins = (x, w) -> bins_equal_number(x, w, num_per_bin=args["num-per-bin"])
    elseif binmethod == "equal-width"
        f_bins = (x, w) -> bins_equal_width(x, w, bin_width=args["bin-width"])
    elseif binmethod == "both"
        f_bins = (x, w) -> bins_both(x, w, num_bins_min=args["num-bins"], bin_width_min=args["bin-width"])
    else
        throw(ArgumentError("bin-method must be one of 'equal-number', 'equal-width', or 'both'. Got $binmethod."))
    end

    if logarithmic == true
        bins = (x, w) -> 10 .^ f_bins(log10.(x), w)
    else
        bins = f_bins
    end

    return bins
end


"""
    add_bin_args(settings)

Returns ArgParseSettings settings with the following binning options added.

- `bin-method`
- `num-per-bin` Effective number of observations in each bin.
- `bin-width` The width of the bin
"""
function add_bin_args(settings::ArgParseSettings)

    hist_settings = ArgParseSettings()

    @add_arg_table hist_settings begin
        "--bin-method"
            help = "the method to use for binning"
            default = "equal-width"
        "--num-per-bin"
            help = "the (minimum) number of points per bin"
            arg_type=Int
        "--bin-width"
            help = "the (minimum) bin width to use"
            arg_type=Float64
    end

    settings = import_settings(settings, hist_settings)

    return settings
end


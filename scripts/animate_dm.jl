#!/usr/bin/env julia
using ArgParse
using LilGuys
using CairoMakie
using Printf
using HDF5
using Arya

# Define default colors
const DEFAULT_COLORS = [
    Makie.RGBA{Float32}(1.0f0, 1.0f0, 1.0f0, 1.0f0),      # White
    COLORS[7], # yellow
    COLORS[2], # orange
    COLORS[6], # rust
    COLORS[4], # pink
]


function get_args()
    s = ArgParseSettings(
         description="Animates the dark matter from a (projected) bird's eye view",
        version="0.2.0"
    )

    @add_arg_table s begin
        "--input", "-i"
            help="Input HDF5 files containing 2D densities"
            required=true
            nargs='+'
        "--colors", "-c"
            help="Colors for each input file (e.g., red,green,blue). Defaults to predefined colors."
            nargs='+'
        "--scalings", "-s"
            help="Scaling factors for each input file. Defaults to 1.0 for each. May be a float, `max` or `match#` where # is the index of the file to match the scaling to."
            nargs='+'
        "-o", "--output"
            help="Output directory for animation frames. Defaults to figures/dm_animation/"
            default="figures/dm_animation/"
        "-P", "--power"
            help="Power for potential scaling. (E.g. 1/2 to scale luminosity as sqrt of density). Note that the power is applied before the intensities are added. Set to 0 for logarithmic"
            default=1
            arg_type=Float64
        "--scalebar"
            help="Length of scalebar in kpc. Setting to 0 disables it"
            default=50
            arg_type=Float64
    end

    args = parse_args(s)
    validate_args!(args)
    return args
end



function validate_args!(args)
    # Validate number of colors and scalings
    num_inputs = length(args["input"])
    
    if !isempty(args["colors"]) && length(args["colors"]) != num_inputs
        error("Number of colors must match number of input files.")
    end

    if !isempty(args["scalings"]) && length(args["scalings"]) != num_inputs
        error("Number of scalings must match number of input files.")
    end

    # Assign default colors if not enough colors provided
    if isempty(args["colors"])
        args["colors"] = DEFAULT_COLORS[1:num_inputs]
    elseif length(args["colors"]) < num_inputs
        error("Not enough colors provided.")
    end

    # Assign default scalings if not enough scalings provided
    if isempty(args["scalings"])
        args["scalings"] = ["max"; fill("match1", num_inputs - 1)...]
        println(args["scalings"])
    end
end

struct Scaling
    factor::Float64
    match::Int
    max::Bool 
end


function Scaling(s::String)
    float_regex = r"^(([0-9]*[.])?[0-9]+)"

    if contains(s, "match")
        index = match(r"match(\d+)", s).captures[1]
        index = parse(Int, string(index))
        factor = match(r"^(([0-9]*[.])?[0-9]+)match", s)
        if factor === nothing
            factor = 1
        else
            factor = parse(Float64, factor.captures[1])
        end
        return Scaling(factor, index, false)
    elseif contains(s, "max")
        factor = match(r"^(([0-9]*[.])?[0-9]+)max", s)
        if factor === nothing
            factor = 1
        else
            factor = parse(Float64, factor.captures[1])
        end
        return Scaling(factor, 0, true)
    else
        return Scaling(parse(Float64, s), 0, false)
    end
end


function parse_scalings(scalings, h5files)
    scalings = Scaling.(scalings)

    scalings_parsed = zeros(Float64, length(scalings))

    for i in 1:length(scalings)
        scaling = scalings[i]
        if scaling.max
            scalings_parsed[i] = scaling.factor / get_max_density(h5files[i])
            @info "scaling $i = max"
        elseif scaling.match > 0
            # pass
        else
            scalings_parsed[i] = scaling.factor
        end
    end

    for i in 1:length(scalings)
        scaling = scalings[i]
        if scaling.match > 0
            scalings_parsed[i] = scaling.factor * scalings_parsed[scaling.match]
        end
        @info "scaling $i = $(scalings_parsed[i])"
    end

    return scalings_parsed
end


function is_static(h5file)
    if "density" ∈ keys(h5file)
        return true
    else
        return false
    end
end


"""
Returns the maximum density value in the given HDF5 file
"""
function get_max_density(h5file)
    if is_static(h5file)
        return maximum(h5file["density"][:, :])
    else
        return maximum([maximum(h5file["/$(key)/density"][:, :]) for key in keys(h5file)])
    end
end


"""
Checks that all files have the same x and y bins.

TODO: eventually allow for ability to just rescale/resample 
densities.
"""
function have_same_bins(h5files, staticfiles=[])
    key1 = keys(h5files[1])[1]
    xbins1 = h5files[1]["/$(key1)/xbins"][:]
    ybins1 = h5files[1]["/$(key1)/ybins"][:]

    for i in 2:length(h5files)
        for j in 1:length(keys(h5files[i]))
            key = keys(h5files[i])[j]
            xbins = h5files[i]["/$(key)/xbins"][:]
            ybins = h5files[i]["/$(key)/ybins"][:]
            if !( xbins1 ≈ xbins && ybins1 ≈ ybins )
                return false
            end
        end
    end

    for i in 1:length(staticfiles)
        xbins = staticfiles[i]["/xbins"][:]
        ybins = staticfiles[i]["/ybins"][:]
        if !( xbins1 ≈ xbins && ybins1 ≈ ybins )
            return false
        end
    end

    return true
end


function main()
    args = get_args()

    input_files = args["input"]
    colors = args["colors"]
    scalings = args["scalings"]
    output_dir = args["output"]
    dm_power = args["power"]
    scalebar_length = args["scalebar"]
    clean_outputs(args["output"])

    # Convert color strings to RGBA
    color_rgba = [Makie.to_color(c) for c in colors]

    # Open all HDF5 files
    h5_files = [h5open(file, "r") for file in input_files]


    try
        static_perm = collect(1:length(h5_files))

        scalings = parse_scalings(scalings, h5_files)

        @info "Checking for static files"
        static_filt = is_static.(h5_files)
        @info "Static files: $(static_filt)"
        static_files = h5_files[static_filt]
        files = h5_files[.!static_filt]
        static_perm = vcat(static_perm[.!static_filt], static_perm[static_filt])
        @info "Static perm: $static_perm"

        @info "checking bins"
        if !have_same_bins(files, static_files)
            error("All input files must have the same bins.")
        end

        scalings = scalings[static_perm]
        color_rgba = color_rgba[static_perm]
        animate_multiple(files, static_files, output_dir;
            scalebar=scalebar_length,
            dm_power=dm_power,
            colors=color_rgba,
            scalings=scalings
        )
    finally
        # Ensure all files are closed
        for f in h5_files
            close(f)
        end
    end
end



function clean_outputs(output_dir::String)
    if isdir(output_dir)
        @info "Removing old output directory: $output_dir"
        rm(output_dir; recursive=true)
    end

    mkpath(output_dir)
end


function animate_multiple(files::Vector{HDF5.File}, static_files::Vector{HDF5.File}, animation_dir::String;
         colors::Vector{Makie.RGBA{Float32}},
         scalings::Vector{Float64},
         scalebar::Float64,
         dm_power::Float64,
    )
    ks_sorted = get_sorted_keys(files)

    Σ_static = [static_files[i]["density"][:, :] * scalings[length(files) + i] for i in eachindex(static_files)]

    for frame in eachindex(ks_sorted)
        print("Creating frame $frame / $(length(ks_sorted))\r")

        xrange, yrange, rgb_image = combine_densities(files, Σ_static, scalings, colors, ks_sorted[frame]; dm_power=dm_power)

        # Create frame
        fig = Figure(figure_padding=(0,0,0,0), backgroundcolor=:black, size=(200,200))
        ax = Axis(fig[1,1], 
            aspect = DataAspect(),
         )

        image!(ax, xrange, yrange, rgb_image)

        hidespines!(ax)
        resize_to_layout!(fig)
        hidedecorations!(ax)

        # Add scalebar if needed
        if scalebar > 0
            add_scalebar!(ax, xrange, yrange, scalebar)
        end

        # Save frame
        Makie.save(joinpath(animation_dir, "frame_$(frame).png"), fig)
    end
end


function combine_densities(files, Σ_static, scalings, colors, key; dm_power=0)
    # Initialize RGB arrays
    combined_R = zeros(Float32, size(files[1]["/$(key)/density"][:,:]))
    combined_G = zeros(Float32, size(files[1]["/$(key)/density"][:,:]))
    combined_B = zeros(Float32, size(files[1]["/$(key)/density"][:,:]))

    local xbins, ybins

    for i in eachindex(files)
        density = files[i]["/$(key)/density"][:, :] * scalings[i]
        # Map intensities with individual color
        color_R, color_G, color_B = map_intensities(density, colors[i])

        xbins = files[i]["/$(key)/xbins"][:]
        ybins = files[i]["/$(key)/ybins"][:]
        # Combine RGB channels
        combined_R .+= color_R
        combined_G .+= color_G
        combined_B .+= color_B
    end

    for i in eachindex(Σ_static)
        density = Σ_static[i]
        # Map intensities with individual color
        color_R, color_G, color_B = map_intensities(density, colors[i + length(files)])
        # Combine RGB channels
        combined_R .+= color_R
        combined_G .+= color_G
        combined_B .+= color_B
    end

    # Normalize combined RGB to [0,1]
    combined_R = clamp.(combined_R, 0.0f0, 1.0f0)
    combined_G = clamp.(combined_G, 0.0f0, 1.0f0)
    combined_B = clamp.(combined_B, 0.0f0, 1.0f0)

    combined_R = remap_intensity(combined_R, dm_power)
    combined_G = remap_intensity(combined_G, dm_power)
    combined_B = remap_intensity(combined_B, dm_power)
    # Create RGB image
    rgb_image = RGBf.(combined_R, combined_G, combined_B)

    return extrema(xbins), extrema(ybins), rgb_image
end


"""
Remaps intensity values to a new range based on a power function.
"""
function remap_intensity(X, power)
    if power < 0
        logmin = power
        logmax = 0

        X = log10.(X)
        X = (max.(X, logmin) .- logmin) ./ (logmax - logmin)
    else
        X = X .^ power
    end

    return X
end


function get_sorted_keys(files)
    num_files = length(files)
    ks_list = [collect(keys(f)) for f in files]

    # Assume all files have the same keys
    for i in 2:num_files
        if length(ks_list[i]) != length(ks_list[1]) || any(ks_list[i] .!= ks_list[1])
            error("All input files must have the same snapshots.")
        end
    end

    ks = ks_list[1]
    idxs = [read_attribute(files[1][key], "snapshot") for key in ks]
    sorted_order = sortperm(idxs)
    ks_sorted = ks[sorted_order]
    idxs_sorted = sort(idxs)
    return ks_sorted
end


function add_scalebar!(ax, xrange, yrange, scalebar::Float64; font="Arial", fontsize=7.5)
    x1 = xrange[1] + scalebar / 2
    y1 = yrange[1] + scalebar / 2
    x2 = x1 + scalebar
    y2 = y1

    lines!(ax, [x1, x2], [y1, y2], color=:grey)

    # Format scalebar label
    label = @sprintf("%.0f kpc", scalebar)
    text!(ax, x1, y1, text=label, color=:grey, font=font, fontsize=fontsize)
end


# Function to map intensities to RGB channels based on a single color
function map_intensities(Sigma_dm, color::Makie.RGBA{Float32})
    color = normalize_color(color)
    color_R = I_to_SI.(clamp.(Sigma_dm .* color.r, 0.0f0, 1.0f0))
    color_G = I_to_SI.(clamp.(Sigma_dm .* color.g, 0.0f0, 1.0f0))
    color_B = I_to_SI.(clamp.(Sigma_dm .* color.b, 0.0f0, 1.0f0))
    return color_R, color_G, color_B
end

# Normalize color to ensure it's within [0,1]
function normalize_color(color::Makie.RGBA{Float32})
    return Makie.RGBA{Float32}(clamp(color.r, 0.0f0, 1.0f0),
                              clamp(color.g, 0.0f0, 1.0f0),
                              clamp(color.b, 0.0f0, 1.0f0),
                              clamp(color.alpha, 0.0f0, 1.0f0))
end

"""
Map luminosity intensity to sRGB scaling
"""
function I_to_SI(I)
    if isnan(I)
        return 0.0f0
    elseif I == Inf
        return 1.0f0
    end

    @assert 0.0f0 <= I <= 1.0f0 "$I is out of range"

    if I < 0.03928 / 12.92
        return I * 12.92
    else
        return I^(1/2.4) * 1.055 - 0.055
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


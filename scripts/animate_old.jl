#!/usr/bin/env julia
using ArgParse

using LilGuys
using CairoMakie

using Printf
using HDF5

const color_0 = Makie.RGBA{Float32}(1.0f0, 1.0f0, 1.0f0, 1.0f0)
const color_1 = Makie.RGBA{Float32}(0.9411765f0,0.89411765f0,0.25882354f0,1.0f0)

function get_args()
    s = ArgParseSettings(
         description="Animates the dark matter from a (projected) bird's eye view",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
            help="Input HDF5 file conaining 2D densities"
            default="."
        "-o", "--output"
            help="Output directory for animation frames. Defaults to figures/dm_animation/"
            default="figures/dm_animation/"
        "-P", "--power"
            help="Power for potential scaling. (E.g. 1/2 to scale luminosity as sqrt of density). Set to 0 for logarithmic"
            default=1
            arg_type=Float64
        "--scalebar"
            help="length of scalebar in kpc. Setting to 0 disables it"
            default=50
            arg_type=Float64
    end

    args = parse_args(s)

    return args
end



function main()
    args = get_args()

    if isdir(args["output"])
        @info "removing old output"
        rm(args["output"]; recursive=true)
    end

    mkpath(args["output"])

    h5open(args["input"], "r") do f
        animate(f, args["output"];
            scalebar=args["scalebar"],
            dm_power=args["power"],
           )
    end
end


function animate(file, animation_dir;
        colorrange=nothing,
        kwargs...
    )
    ks = collect(keys(file))


    idxs = [read_attribute(file[key], "snapshot") for key in ks]
    ks = ks[sortperm(idxs)]
    idxs = sort(idxs)

    if colorrange == nothing
        @info "finding maximum color"
        colormax = maximum([maximum(file["/$(key)/density"][:,:]) for key in ks])
        colorrange = (0, colormax)
    end

    @info "using colorrange: $colorrange"

    for frame in eachindex(idxs)
        print("creating frame $frame / $(length(idxs))\r")

        density = file["/$(ks[frame])/density"][:, :]
        xbins= file["/$(ks[frame])/xbins"][:]
        ybins = file["/$(ks[frame])/ybins"][:]

        fig, ax, p = make_frame(xbins, ybins, density;
            colorrange=colorrange, 
            kwargs...
        )
        Makie.save(joinpath(animation_dir, "frame_$frame.png"), fig)
    end
    println()
end


function make_frame(xbins, ybins, density; colorrange, font="Arial", scalebar=50, Sigma_mw=nothing,
        color_A=color_0, color_B=color_1,
        dm_power = 1/5
    )
    fig = Figure(figure_padding=(0,0,0,0), size=(200,200), backgroundcolor=:black) # px scale is x5 so this works for 1000 bins
    ax = Axis(fig[1, 1],
        aspect=DataAspect(),
    )
    
    xrange, yrange = extrema(xbins), extrema(ybins)

    if dm_power == 0
        Sigma_dm = log10.(density./ colorrange[2])
    else
        Sigma_dm = (density./ colorrange[2]) .^ dm_power
    end

    colors = map_intensities(Sigma_dm, color_A)

    p = image!(xrange, yrange, colors)

    hidespines!(ax)
    resize_to_layout!(fig)
    hidedecorations!(ax)

    # create scalebar

    if scalebar > 0
        x1 = xrange[1] + scalebar / 2
        y1 = yrange[1] + scalebar / 2
        x2 = x1 + scalebar
        y2 = y1

        lines!([x1, x2], [y1, y2], color=:grey)

        # if reasonable to round, than do
        if (scalebar % 1 < 0.001) || (scalebar % 1 > 0.999)
            scalebar = @sprintf("%.0f", scalebar)
        elseif (scalebar % 0.1 < 0.0001) || (scalebar % 0.1 > 0.0999)
            scalebar = @sprintf("%.1f", scalebar)
        else
            scalebar = @sprintf("%.2f", scalebar)
        end
        text!(x1, y1, text="$scalebar kpc", color=:grey, font=font, fontsize=7.5)
    end

    return Makie.FigureAxisPlot(fig, ax, p)
end



function I_to_SI(I)
    if isnan(I)
        return 0.
    end
    if I === Inf
        return 1.
    end
    
    @assert 0 <= I <= 1 "$I is out of range"

    if I < 0.03928 / 12.92
        return I * 12.92
    else
        return I^(1/2.4) * 1.055 - 0.055
    end
end

function normalize_colors(colors)
    if isnan(colors)
        return 0.0
    elseif colors === Inf
        return 1.0
    end
    colors = min(colors, 1.0)
    colors = max(colors, 0.0)
    
    return colors
end

function map_intensities(A, B, color_A, color_B)
    R = @. A * Makie.red(color_A) + B * Makie.red(color_B)
    G = @. A * Makie.green(color_A) + B * Makie.green(color_B)
    B = @. A * Makie.blue(color_A) + B * Makie.blue(color_B)
    R = normalize_colors.(R) 
    G = normalize_colors.(G) 
    B = normalize_colors.(B)
    color_R = @. I_to_SI.(R)
    color_G = @. I_to_SI.(G)
    color_B = @. I_to_SI.(B)

    
    return RGBf.(color_R, color_G, color_B)
end

function map_intensities(A, color_A;)
    color_R = I_to_SI.(normalize_colors.(A .* Makie.red.(color_A)))
    color_G = I_to_SI.(normalize_colors.(A .* Makie.green.(color_A)))
    color_B = I_to_SI.(normalize_colors.(A .* Makie.blue.(color_A)))

    
    return RGBf.(color_R, color_G, color_B)
end







if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

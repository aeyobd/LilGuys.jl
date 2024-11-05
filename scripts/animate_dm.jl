#!/usr/bin/env julia
using ArgParse

using LilGuys
using CairoMakie

using StatsBase
using PythonCall
agama = pyimport("agama")
using Printf


# TODO: could add depth by being careful with the MW transparency ordering...


const color_0 = Makie.RGBA{Float32}(0.9411765f0,0.89411765f0,0.25882354f0,1.0f0)

function get_args()
    s = ArgParseSettings(
         description="Animates the dark matter from a (projected) bird's eye view",
        version="0.1.0"
    )

    @add_arg_table s begin
        "input"
            help="Input (simulation output) file"
            default="."

        "-o", "--output"
            help="Output directory for animation frames. Defaults to figures/dm_animation/"
            default="figures/dm_animation/"
        "-k", "--skip"
            help="Skip length for snapshots"
            default=10
            arg_type=Int
        "-t", "--test"
            help="Test mode: only process the first and last snapshot"
            action="store_true"
        "--limits"
            help="limits for plot"
            nargs='+'
            arg_type=Float64
            default=[200]
        "--x_vec"
            help="x vector for projection"
            nargs=3
            arg_type=Float64
            default=[0, 1, 0]
        "--y_vec"
            help="y vector for projection"
            nargs=3
            arg_type=Float64
            default=[0, 0, 1]
        "--n_bins"
            help="number of bins for histogram"
            default=1000
            arg_type=Int
        "-p", "--potential"
            help="stellar potential file for MW representation. If set to nothing, then do not add disk potential to animation"
        "--font"
            help="font for scalebar"
            default="/astro/dboyea/fonts/Arev.ttf"
    end

    args = parse_args(s)

    if length(args["limits"]) == 1
        args["limits"] = (-args["limits"][1], args["limits"][1], -args["limits"][1], args["limits"][1])
    elseif length(args["limits"]) == 4
        args["limits"] = (args["limits"][1], args["limits"][2], args["limits"][3], args["limits"][4])
    else
        @error "May only specify one or four limits"
    end

    return args
end




function main()
    args = get_args()

    if isdir(args["output"])
        rm(args["output"]; recursive=true)
    end

    mkpath(args["output"])

    @info "loading sample"
    out = Output(args["input"])
    
    xbins = LinRange(args["limits"][1], args["limits"][2], args["n_bins"])
    ybins = LinRange(args["limits"][3], args["limits"][4], args["n_bins"])
    bins = (xbins, ybins)
    @assert issorted(xbins) && issorted(ybins)

    if args["test"]
        idx = [1, length(out)]
    else
        idx = 1:args["skip"]:length(out)
    end

    animate(out, bins, args["output"]; 
            idx=idx, 
            x_vec=args["x_vec"], y_vec=args["y_vec"], 
            Sigma_disk=args["potential"],
            font=args["font"]
           )

end


function get_xy(out, idx; x_vec=[0, 1, 0], y_vec = [0, 0, 1])
    # shortcut for hdf5
    # -1 since HDF5 is zero indexed
    idx -= 1
    pos = out.h5file["snap$idx/PartType1/Coordinates"][:, :]
    A = [x_vec y_vec]'
    xy = A * pos
    return xy[1, :], xy[2, :]
end


function project_points(x, y, xybins)
    h1 = fit(Histogram, (x, y), xybins )
end

function transparency_map(x) 
    x = x
    Makie.RGBAf(color_0.r, color_0.g, color_0.b * x, x)
end



function make_frame(h; colorrange, font, scalebar=50, Sigma_disk=nothing)
    fig = Figure(figure_padding=(0,0,0,0), size=(200,200)) # px scale is x5 so this works for 1000 bins
    ax = Axis(fig[1, 1],
        aspect=DataAspect(),
    )
    

    xrange, yrange = extrema.(h.edges)

    p = image!(xrange, yrange, log10.(h.weights), colorrange=colorrange)

    if Sigma_disk !== nothing
        image!(xrange, yrange, transparency_map.(normalized_density(Sigma_disk)))
    end
    hidespines!(ax)
    resize_to_layout!(fig)
    hidedecorations!(ax)

    # create scalebar

    x1 = xrange[1] + scalebar / 2
    y1 = yrange[1] + scalebar / 2
    x2 = x1 + scalebar
    y2 = y1

    lines!([x1, x2], [y1, y2], color=:grey)
    text!(x1, y1, text="$scalebar kpc", color=:grey, font=font, fontsize=7.5)

    return Makie.FigureAxisPlot(fig, ax, p)
end


function animate(out, bins, animation_dir; 
        idx=eachindex(out), 
        x_vec=[0, 1, 0], y_vec = [0, 0, 1], colorrange=nothing,
        Sigma_disk=nothing, scalebar=50, font="Arial"
    )

    if colorrange == nothing
        x, y = get_xy(out, 1; x_vec=x_vec, y_vec=y_vec)
        h = project_points(x, y, bins)
        colorrange = [0, maximum(log10.(h.weights))]
    end

    @info "using colorrange: $colorrange"

    for (n, i) in enumerate(idx)
        print("creating frame $n / $(length(idx))\r")
        x, y = get_xy(out, i; x_vec=x_vec, y_vec=y_vec)
        h = project_points(x, y, bins)
        fig, ax, p = make_frame(h, colorrange=colorrange, scalebar=scalebar,
            font=font
        )
        Makie.save(joinpath(animation_dir, "frame_$i.png"), fig)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

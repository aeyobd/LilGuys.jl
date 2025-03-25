module MakieExt

using Makie
using StatsBase

import Makie: convert_arguments, convert_single_argument
using LilGuys
import LilGuys: @assert_3vector
import LinearAlgebra: norm, dot

using Arya

# methods defined here
import LilGuys: plot_xyz, plot_xyz!
import LilGuys: cmd_axis
import LilGuys: projecteddensity, projecteddensity!
import LilGuys: hide_grid!, plot_log_Σ!, plot_Γ!
import LilGuys: @savefig



const log_R_label = L"\log\,(r\, /\, \mathrm{kpc})"
const log_rho_label = L"\log\,(\rho\, /\, 10^{10}\mathrm{M}_\odot\mathrm{pc}^{-3})"
const v_circ_label = L"{v}_\mathrm{circ}\, /\, \mathrm{km\, s}^{-1}"
const xi_arcmin_label = L"\xi\, /\, \mathrm{arcmin}"
const eta_arcmin_label = L"\eta\, /\, \mathrm{arcmin}"

const markers = [:circle, :rect, :star5, :cross, :diamond, :hexagon, :octagon, :star4, :star6, :star7, :star8, :star9, :star10, :star11, :star12, :star13, :star14, :star15, :star16, :star17, :star18, :star19, :star20]



"""
    projecteddensity(snapshot)

Projects a snapshot into 2 dimensions and plots the density.

# Attributes
- `snapshot`: The snapshot to plot.
- `bin=200`: The number of bins to use.
- `r_max=10`: The maximum radius to plot.
- `centre=true`: If true, centres the plot on the centre of mass. Otherwise, plots in the range -r_max to r_max.
- `xdirection=1`: The direction to plot in the x axis.
- `ydirection=2`: The direction to plot in the y axis.

"""
@recipe(ProjectedDensity, snapshot) do scene
    Attributes(
        bins = 200,
        r_max = 10,
        centre = true,
        xdirection = 1,
        ydirection = 2,
        colorrange = theme(scene, :colorrange),
        colormap = theme(scene, :colormap),
        colorscale = identity,
    )
end

Makie.needs_tight_limits(::ProjectedDensity) = true

function Makie.plot!(p::ProjectedDensity)
    snap = p[:snapshot][]
    
    # Extract positions and calculate limits
    x0 = snap.x_cen[p[:xdirection][]]
    y0 = snap.x_cen[p[:ydirection][]]
    r_max = p[:r_max][]
    
    limits = p[:centre][] ? 
        (-r_max + x0, r_max + x0, -r_max + y0, r_max + y0) : 
        (-r_max, r_max, -r_max, r_max)
    
    # Extract positions and plot
    x = snap.positions[p[:xdirection][], :]
    y = snap.positions[p[:ydirection][], :]
    
    # Use Arya's hist2d! function
    Arya.hist2d!(p, x, y; 
        bins=p[:bins][], 
        limits=limits, 
        weights=snap.masses,
        colorrange=p[:colorrange][],
        colormap=p[:colormap][],
        colorscale=p[:colorscale][]
    )
    
    return p
end


"""
    plot_xyz(args...; plot!, labels, units, limits, kwargs...)


Given (any number of) 3xN matricies of xyz positions, makes orbit plots in each plane.

# Arguments
- `args...`: 3xN matricies of xyz positions.
- `plot`: The plotting mode to use. May be :scatter or :lines
- `labels`: Labels for each orbit.
- `idx_scatter`: If set, also plots points at the given index for each 
- `xlabel`, `ylabel`, `zlabel`: Labels for the axes.
- `units`: Units for the axes.
- `limits`: Limits for the axes. Should be a tuple of tuples for x y and z limits
- `kwargs...`: Additional keyword arguments to pass to the plotting function.
"""
function plot_xyz(args...; 
        labels=nothing, 
        limits=nothing, 
        idx_scatter=nothing, 
        times=nothing, 
        xlabel="x", ylabel="y", zlabel="z", 
        units=" / kpc", 
        colorrange=nothing, 
        kwargs...
    )

    limits = limits_xyz(args...; limits=limits)
    fig = Figure()
    axes = axis_xyz(fig; limits=limits, xlabel=xlabel, ylabel=ylabel, 
                    zlabel=zlabel, units=units)

    if colorrange === nothing && times !== nothing
        colorrange = extrema(times * T2GYR)
        kwargs = (; colorrange=colorrange, kwargs...)
    end

    p = plot_xyz!(axes, args...; labels=labels, times=times, kwargs...)

    if idx_scatter !== nothing
        args_scatter = [arg[:, idx] for (arg, idx) in zip(args, idx_scatter)]

        if times !== nothing
            times = times[idx_scatter[1]]
        else
            times = nothing
        end 


        plot_xyz!(axes, args_scatter...; times=times, plot=:scatter, labels=nothing, kwargs...)
    end

    if labels !== nothing
        Legend(fig[1, 2], axes[1], tellwidth=false)
    elseif times !== nothing
        cbar = Colorbar(fig[1, 2], p, label="time / Gyr", tellwidth=false, halign=:left)
    end


    return fig
end



function plot_xyz!(axes, args...; plot = :lines, labels=nothing, 
        times=nothing, limits=nothing, kwargs...)

    if plot === :lines
        plot! = lines!
    elseif plot === :scatter
        plot! = scatter!
    else
        error("Unknown plot type: $plot")
    end

    ax_xy, ax_yz, ax_xz = axes

    Nargs = length(args)

    if times !== nothing
        kwargs = (; color=times * T2GYR, kwargs...)
    end

    for i in 1:Nargs
        @assert_3vector args[i]
    end

    local p

    for i in 1:Nargs
        if labels !== nothing
            label = labels[i]
        else 
            label = nothing
        end

        plot!(ax_xy, args[i][1, :], args[i][2, :]; label=label, kwargs...)
        plot!(ax_yz, args[i][2, :], args[i][3, :]; label=label, kwargs...)
        p = plot!(ax_xz, args[i][1, :], args[i][3, :]; label=label, kwargs...)
    end



    return p
end


"""
    axis_xyz(fig; limits, xlabel="x", ylabel="y", zlabel="z", units=" / kpc")

Creates 3 linked axes for a coupled xyz plot of a 3D system.
"""
function axis_xyz(fig; limits, 
        xlabel="x", ylabel="y", zlabel="z", units=" / kpc")
    ax_xy = Axis(fig[1, 1], xlabel="$xlabel$units", ylabel="$ylabel$units", 
                 aspect=DataAspect(), limits=limits[1:2])
    ax_yz = Axis(fig[2, 2], xlabel="$ylabel$units", ylabel="$zlabel$units", 
                 aspect=DataAspect(), limits=limits[2:3])
    ax_xz = Axis(fig[2, 1], xlabel="$xlabel$units", ylabel="$zlabel$units", 
                 aspect=DataAspect(), limits=limits[[1, 3]])

    # hide extra axis labels since they are linked
    linkxaxes!(ax_xy, ax_xz)
    hidexdecorations!(ax_xy, grid=false, ticks=false, minorticks=false)
    linkyaxes!(ax_xz, ax_yz)
    hideydecorations!(ax_yz, grid=false, ticks=false, minorticks=false)

    # fixes spacing...
    rowsize!(fig.layout, 1, Aspect(1, 1))
    rowsize!(fig.layout, 2, Aspect(1, 1))

    resize_to_layout!(fig)

    return ax_xy, ax_yz, ax_xz
end



"""
    limits_xyz(args...; symmetric=true, limits=nothing)

Calculates limits for a 3D plot based on the maximum and minimum values of the input arrays.
If `symmetric=true`, the limits are symmetric about the origin. Otherwise, the limits are set to the maximum extent of the data and centred on the centre of the data range.
"""
function limits_xyz(args...; symmetric=true, limits=nothing)
    for i in eachindex(args)
        @assert_3vector args[i]
    end

    if limits === nothing
        xmax = maximum(maximum.(args))
        xmin = minimum(minimum.(args))

        ymax = maximum(maximum.(args))
        ymin = minimum(minimum.(args))

        zmax = maximum(maximum.(args))
        zmin = minimum(minimum.(args))

        if symmetric
            xmax = max(abs(xmax), abs(xmin))
            ymax = max(abs(ymax), abs(ymin))
            zmax = max(abs(zmax), abs(zmin))

            rmax = 1.1 * max(xmax, ymax, zmax)

            limits = ((-rmax, rmax), (-rmax, rmax), (-rmax, rmax))
        else
            rmax = 1/2 * max(xmax - xmin, ymax - ymin, zmax - zmin)
            xm = (xmax + xmin) / 2
            ym = (ymax + ymin) / 2
            zm = (zmax + zmin) / 2

            limits = ((xm - rmax, xm + rmax), (ym - rmax, ym + rmax), (zm - rmax, zm + rmax))
        end
    end

    return limits
end





function convert_arguments(p::Type{<:Scatter}, h::LilGuys.StellarDensityProfile)
    return (h.log_R, LilGuys.middle.(h.log_Sigma))
end


function convert_single_argument(a::AbstractArray{T}) where {T<:LilGuys.Measurement}
    return LilGuys.middle.(a)
end


"""
    plot_log_Σ!(ax, p; kwargs...)

Plot a density profile from a LilGuys.StellarDensityProfile.
kwargs passed to Arya.errorscatter!
"""
function plot_log_Σ!(ax, p::LilGuys.StellarDensityProfile; kwargs...)
    x = p.log_R
    y = LilGuys.middle.(p.log_Sigma)
    yerror = LilGuys.credible_interval.(p.log_Sigma)
    filt = isfinite.(y)

    x = x[filt]
    y = y[filt]
    yerror = yerror[filt]

    errorscatter!(ax, x, y; yerror=yerror,  kwargs...)
end



"""
    plot_Γ!(ax, p; kwargs...)

Plot a density profile from a LilGuys.StellarDensityProfile.
kwargs passed to Arya.errorscatter!
"""
function plot_Γ!(ax, p::LilGuys.StellarDensityProfile; log_R=true, kwargs...)
    x = p.log_R
    y = LilGuys.middle.(p.Gamma)
    yerror = LilGuys.credible_interval.(p.Gamma)
    filt = isfinite.(y)

    x = x[filt]
    y = y[filt]
    yerror = yerror[filt]

    if !(log_R)
        x = 10 .^ x
    end

    errorscatter!(ax, x, y; yerror=yerror,  kwargs...)
end



"""
    hide_grid!(ax)

Hides the grid lines on an axis.
"""
function hide_grid!(ax)
    ax.xgridvisible = false
    ax.ygridvisible = false
end


"""
    @savefig name [fig=fig]

Saves a figure to the current figdir as both a pdf and a jpeg.
If `fig` is not provided, it defaults to the current figure (assumed to be `fig`).
`FIGDIR` may be defined as the path where figures are saved.
`FIGSUFFIX` may be defined as the suffix for figures (e.g. notebook name).
Saves figures with the basename name and in both pdf and png formats.
Returns the figure object.

# Notes
Will make the FIGDIR if this does not exist (but will not mkpath).
"""
macro savefig(name, fig=nothing)
    if isnothing(fig)
        fig = esc(:fig)
    else
        fig = esc(fig)
    end


    return quote
        local dir, filepath, filename

        if isdefined(@__MODULE__, :FIGDIR)
            dir = (@__MODULE__).FIGDIR
            if !isdir(dir)
                @info "Creating figure directory $dir"
                mkdir(dir)
            end
        else
            dir = ""
        end

        filename = $(esc(name))
        if isdefined(@__MODULE__, :FIGSUFFIX)
            suffix = (@__MODULE__).FIGSUFFIX
            filename *= "$(suffix)"
        end

        filepath = joinpath(dir, filename)


        @info "Saving figure to $filename.pdf and ---.png in $dir"
        Makie.save(filepath * ".pdf", $fig, pt_per_unit=1)
        Makie.save(filepath * ".png", $fig, pt_per_unit=1)

        $fig
    end 
end


end # module 

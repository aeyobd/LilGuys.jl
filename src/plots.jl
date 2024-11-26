module Plots

using Makie
using StatsBase

using ..LilGuys
using Arya

import LinearAlgebra: norm, dot


export hide_grid!

const log_r_label = L"\log\,(r\, /\, \mathrm{kpc})"
const log_rho_label = L"\log\,(\rho\, /\, 10^{10}\mathrm{M}_\odot\mathrm{pc}^{-3})"
const v_circ_label = L"{v}_\mathrm{circ}\, /\, \mathrm{km\, s}^{-1}"
const xi_arcmin_label = L"\xi\, /\, \mathrm{arcmin}"
const eta_arcmin_label = L"\eta\, /\, \mathrm{arcmin}"

const markers = [:circle, :rect, :star5, :cross, :diamond, :hexagon, :octagon, :star4, :star6, :star7, :star8, :star9, :star10, :star11, :star12, :star13, :star14, :star15, :star16, :star17, :star18, :star19, :star20]


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


function limits_xyz(args...; limits=nothing)
    Nargs = length(args)
    if limits === nothing
        xmax = maximum([maximum(abs.(args[i][1, :])) for i in 1:Nargs])
        ymax = maximum([maximum(abs.(args[i][2, :])) for i in 1:Nargs])
        zmax = maximum([maximum(abs.(args[i][3, :])) for i in 1:Nargs])
        # a little bit of breathing room.
        rmax = 1.1 * max(xmax, ymax, zmax)
        limits = ((-rmax, rmax), (-rmax, rmax), (-rmax, rmax))
    end

    return limits
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
function plot_xyz(args...; labels=nothing, limits=nothing, idx_scatter=nothing, times=nothing,
        xlabel="x", ylabel="y", zlabel="z", units=" / kpc", colorrange=nothing, kwargs...)

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
        println(args_scatter)


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
    plot_centre!(positions; rotation, width, z_filter, marker_z, label, kwargs...)

Given a 3xN matrix of xyz positions, plots only points beloging to the very centre.
"""
function plot_centre!(positions; rotation=(0,0), width=5, z_filter=:sphere, 
        marker_z=nothing, label="", kwargs...)
    if z_filter === :sphere
        filt = calc_r(positions) .< width
    elseif z_filter === :box
        filt = abs.(positions[3, :]) .< width
    elseif z_filter === :none
        filt = trues(size(positions, 2))
    else
        error("Unknown z_filter: $z_filter")
    end

    x = positions[1, filt]
    y = positions[2, filt]
    z = positions[3, filt]

    if marker_z === "z"
        marker_z = z
    elseif marker_z !== nothing
        marker_z = marker_z[filt]
    end

    scatter!(x, y; color=marker_z, label=label, kwargs...)
end





"""
Plot the projected 2d density of a snapshot as a heatmap.
"""
function projected_density!(snap; bins=200, r_max=10, centre=true, xdirection=1, ydirection=2, kwargs...)


    x0 = snap.x_cen[1]
    y0 = snap.x_cen[2]
    if centre
        limits = (-r_max + x0, r_max + x0, -r_max + y0, r_max + y0)
    else
        limits = (-r_max, r_max, -r_max, r_max)
    end


    x = snap.positions[xdirection, :]
    y = snap.positions[ydirection, :]

	p = Arya.hist2d!(x, y;
        bins=bins, limits=limits, weights=snap.masses, kwargs...)

    return p
end



function vx_hist_fit!(snap; stars=true, 
        bins=30,
        direction=1,
        limits=nothing,
        kwargs...
    )


    v = snap.velocities[direction, :] * V2KMS
	w = snap.weights

	σ = std(v, weights(w))
    μ = mean(v, weights(w))
    println("μ = $μ, σ = $σ")

    h = Arya.histogram(v, bins, weights=w, normalization=:pdf, limits=limits)
    p = barplot!(h; kwargs...)

    v_model = range(minimum(v), maximum(v), length=100)
    hist_model = LilGuys.gaussian.(v_model, μ, σ)

    lines!(v_model, hist_model; color=:red)


    return p
end


# axis constructors

"""
    cmd_axis(gs; kwargs...)

A gaia CMD axis
"""
function cmd_axis(gs::Makie.GridPosition, kwargs...)
    ax = Axis(gs; 
              xlabel="Bp - Rp (mag)", 
              ylabel="G (mag)", 
              yreversed=true, 
              kwargs...
        )

    return ax
end


"""
    xy_axis(gs; xlabel, ylabel, units, kwargs...)

A generic axis for plotting x-y coordinate (orthoganaly projected) data.
"""
function xy_axis(gs; xlabel="x", ylabel="y", units="specify", kwargs...)
    ax = Axis(gs;
        aspect=DataAspect(), 
        xlabel = "$xlabel / $units", 
        ylabel="$ylabel / $units",
        kwargs...
    )

    return fig, ax
end



"""
    rho_axis(gs; kwargs...)

The axis for plotting the density profile.
"""
function rho_axis(gs; kwargs...)
    ax = Axis(gs; 
        xlabel=log_r_label, 
        ylabel=log_rho_label, 
        kwargs...
    )

    return ax
end



"""
    xi_eta_axis(gs; kwargs...)

The tangent plane axis 
"""
function xi_eta_axis(gs; kwargs...)
    ax = Axis(gs; 
        xlabel=xi_arcmin_label, 
        ylabel=eta_arcmin_label,
        aspect=DataAspect(), 
        kwargs...
    )

    return ax
end



"""
    hide_grid!(ax)

Hides the grid lines on an axis.
"""
function hide_grid!(ax)
    ax.xgridvisible = false
    ax.ygridvisible = false
end


end # module plots

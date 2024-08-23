module Plots

using Makie
using StatsBase

using ..LilGuys
using Arya

import LinearAlgebra: norm, dot


const log_r_label = L"\log\,(r\, /\, \mathrm{kpc})"
const log_rho_label = L"\log\,(\rho\, /\, 10^{10}\mathrm{M}_\odot\mathrm{pc}^{-3})"

const v_circ_label = L"{v}_\mathrm{circ}\, /\, \mathrm{km\, s}^{-1}"

"""
Given (any number of) 3xN matricies of xyz positions, makes orbit plots in each plane.
"""
function plot_xyz(args...; plot! =lines!, labels=nothing, units=" / kpc", kwargs...)

    fig = Figure()
    Nargs = length(args)


    ax_xy = Axis(fig[1, 1], xlabel="x$units", ylabel="y$units", aspect=1)
    ax_yz = Axis(fig[2, 2], xlabel="y$units", ylabel="z$units", aspect=1)
    ax_xz = Axis(fig[2, 1], xlabel="x$units", ylabel="z$units", aspect=1)

    for i in 1:Nargs
        if labels !== nothing
            label = labels[i]
        else 
            label = nothing
        end
        plot!(ax_xy, args[i][1, :], args[i][2, :]; label=label, kwargs...)
        plot!(ax_yz, args[i][2, :], args[i][3, :]; label=label, kwargs...)
        plot!(ax_xz, args[i][1, :], args[i][3, :]; label=label, kwargs...)
    end

    linkxaxes!(ax_xy, ax_xz)
    hidexdecorations!(ax_xy, grid=false, ticks=false, minorticks=false)
    linkyaxes!(ax_xz, ax_yz)
    hideydecorations!(ax_yz, grid=false, ticks=false, minorticks=false)

    if labels !== nothing
        Legend(fig[1, 2], ax_xy, tellwidth=false)
    end

    return fig
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




"""
Plot the dark matter density profile of a snapshot.
"""
function plot_ρ_dm!(snap; filt_bound=true, bins=200, kwargs...)
	rs = LilGuys.calc_r(snap)
    mass = snap.masses

    if filt_bound
        filt = LilGuys.get_bound(snap)
        rs = rs[filt]
        mass = mass[filt]
	end

    r, ρ = LilGuys.calc_ρ_hist(rs, bins, weights=mass)
	x = log10.(midpoints(r))
	lines!(x, log10.(ρ); kwargs...)
end



function plot_ρ!(prof::LilGuys.AbstractProfile, rrange=(-2, 2); N=1000, kwargs...)
    log_r = LinRange(rrange[1], rrange[2], N)
    r = 10 .^ log_r
    ρ = LilGuys.calc_ρ.(prof, r)

    lines!(log_r, log10.(ρ); kwargs...)
end


function plot_ρ_s!(snap; bins=200, kwargs...)
	rs = LilGuys.calc_r(snap)
	ps = snap.weights
	r, ρ = LilGuys.calc_ρ_hist(rs, bins, weights=ps)
	x = log10.(midpoints(r))
	lines!(x, log10.(ρ); kwargs...)
end


"""
Plot the circular velocity profile of a snapshot.
"""
function plot_v_circ!(snap; kwargs...)
	rc, vc = LilGuys.calc_v_circ(snap)
	lines!(log10.(rc), V2KMS * vc; kwargs...)
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

function cmd_axis(gs::Makie.GridPosition, kwargs...)
    ax = Axis(gs; xlabel="Bp-Rp", ylabel="G", 
        yreversed=true, kwargs...)
    return ax
end

function xy_axis(directionx=1, directiony=2; kwargs...)
    fig = Figure()

    xlabel,ylabel = ["x", "y", "z"][[directionx, directiony]]

    ax = Axis(fig[1,1], aspect=DataAspect(), 
        xlabel = "$xlabel / kpc", ylabel="$ylabel / kpc"
    )

    return fig, ax
end


function rho_axis(gs; kwargs...)
    ax = Axis(gs; 
        xlabel=log_r_label, ylabel=log_rho_label, kwargs...)
    return ax
end




end # module plots

using Makie



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






using Makie



"""
Given (any number of) 3xN matricies of xyz positions, makes orbit plots in each plane.
"""
function plot_xyz(args...; labels=nothing, units=" / kpc", aspect_ratio=1, kwargs...)

    plots = []

    axis_labels = ["x", "y", "z"]
    Nargs = length(args)

    for (i, j) in [(1, 2), (2, 3), (1, 3)]
        p = plt.plot()
        for k in 1:Nargs
            arg = args[k]
            if labels !== nothing
                label = labels[k]
            else
                label = ""
            end
            plt.plot!(p, arg[i, :], arg[j, :]; aspect_ratio=aspect_ratio, label=label, kwargs...)
        end
        plt.xlabel!(p, "$(axis_labels[i])$units")
        plt.ylabel!(p, "$(axis_labels[j])$units")
        push!(plots, p)
    end

    return plots
end


function plot_xyz_layout(args...; size=(800, 800), labels=nothing, kwargs...)
    p_xy, p_yz, p_xz = plot_xyz(args...; legend=false, labels=labels, kwargs...)

    p_legend = plt.plot(axis=false)
    if labels !== nothing
        for label in labels
            plt.plot!(p_legend, [], [], label=label, legend_position=:topleft)
        end
    end

    layout = plt.@layout([° °; ° ° ])

    return plt.plot(p_xy, p_legend, p_xz, p_yz, layout=layout, size=size)
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

    plt.scatter!(x, y; marker_z=marker_z, label=label, kwargs...)
end


function plot_centre(positions; width=5, kwargs...)
    plt.plot(xlabel="x/kpc", ylabel="y/kpc", aspect_ratio=1,
             xlims=(-width, width), ylims=(-width, width))

    plot_centre!(positions; width, kwargs...)
end

"""
Given (any number of) 3xN matricies of xyz positions, makes orbit plots in each plane.
"""
function scatter_xyz(args...; kwargs...)

    plots = []

    labels = ["x", "y", "z"]

    for (i, j) in [(1, 2), (1, 3), (2, 3)]
        p = plt.plot()
        for arg in args
            plt.scatter!(p, arg[i, :], arg[j, :]; kwargs...)
        end
        plt.xlabel!(p, "$(labels[i]) / kpc")
        plt.ylabel!(p, "$(labels[j]) / kpc")
        push!(plots, p)
    end

    return plots
end





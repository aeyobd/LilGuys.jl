import Plots as plt




"""
Given (any number of) 3xN matricies of xyz positions, makes orbit plots in each plane.
"""
function plot_xyz(args...; labels=nothing, kwargs...)

    plots = []

    axis_labels = ["x", "y", "z"]
    Nargs = length(args)

    for (i, j) in [(1, 2), (1, 3), (2, 3)]
        p = plt.plot()
        for k in 1:Nargs
            arg = args[k]
            if labels !== nothing
                label = labels[k]
            else
                label = ""
            end
            plt.plot!(p, arg[i, :], arg[j, :]; label=label, kwargs...)
        end
        plt.xlabel!(p, "$(axis_labels[i]) / kpc")
        plt.ylabel!(p, "$(axis_labels[j]) / kpc")
        push!(plots, p)
    end

    return plots
end



function plot_centre!(positions; rotation=(0,0), width=5, marker_z=nothing, kwargs...)
    filt = calc_r(positions) .< width
    x = positions[1, filt]
    y = positions[2, filt]
    z = positions[3, filt]

    if marker_z === "z"
        marker_z = z
    elseif marker_z !== nothing
        marker_z = marker_z[filt]
    end

    plt.scatter!(x, y; marker_z=marker_z, kwargs...)
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





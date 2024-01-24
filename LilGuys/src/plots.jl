import Plots as plt




"""
Given (any number of) 3xN matricies of xyz positions, makes orbit plots in each plane.
"""
function plot_xyz(args...; kwargs...)

    plots = []

    labels = ["x", "y", "z"]

    for (i, j) in [(1, 2), (1, 3), (2, 3)]
        p = plt.plot()
        for arg in args
            plt.plot!(p, arg[i, :], arg[j, :]; kwargs...)
        end
        plt.xlabel!(p, "$(labels[i]) / kpc")
        plt.ylabel!(p, "$(labels[j]) / kpc")
        push!(plots, p)
    end

    return plots
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

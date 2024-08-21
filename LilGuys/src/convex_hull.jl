import Polyhedra

"""
    convex_hull(x, y)
Given a vector of x and y coordinates, returns
the convex hull bounding the points.
"""
function convex_hull(xi, eta)
    @assert length(xi) == length(eta)
    filt = .!isnan.(xi) .& .!isnan.(eta)
    @assert sum(filt) > 2

    ps = [[x, e] for (x, e) in zip(xi[filt], eta[filt])]
    p = Polyhedra.convexhull(ps...)
    b = Polyhedra.planar_hull(p).points.points
    return first.(b), last.(b)
end

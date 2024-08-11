import Polyhedra

"""
Given a vector of x and y coordinates, returns
the convex hull bounding the points.
"""
function convex_hull(xi, eta)
    ps = [[x, e] for (x, e) in zip(xi, eta)]
    p = Polyhedra.convexhull(ps...)
    b = Polyhedra.planar_hull(p).points.points
    return first.(b), last.(b)
end

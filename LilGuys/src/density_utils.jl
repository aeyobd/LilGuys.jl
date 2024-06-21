import Base: @kwdef

import LinearAlgebra: diag, dot, norm, normalize

import LilGuys as lguys

using Measurements
import Arya: histogram
import StatsBase: weights, mean
import TOML
import Polyhedra




value(x::Measurement) = x.val
err(x::Measurement) = x.err
value(x) = x
err(x) = zero(x)


F = Float64


"""
An observed 2D density profile
"""
@kwdef mutable struct ObsProfile
    log_r::Vector{F}
    log_r_bins::Vector{F}
    r_units::String
    counts::Vector{F}

    mass_in_annulus::Vector{F}
    mass_in_annulus_err::Vector{F}

    M_in::Vector{F}
    M_in_err::Vector{F}

    Sigma::Vector{F}
    Sigma_err::Vector{F}
    Sigma_m::Vector{F}
    Sigma_m_err::Vector{F}
    log_Sigma::Vector{F}
    log_Sigma_err::Vector{F}

    Gamma::Vector{F}
    Gamma_err::Vector{F}
    Gamma_max::Vector{F}
    Gamma_max_err::Vector{F}
end




"""
general method to convert a struct to a dictionary
"""
function struct_to_dict(S)
    return Dict(key=>getfield(S, key) for key in fieldnames(typeof(S)))
end


function dict_to_tuple(D)
    return NamedTuple((Symbol(key), value) for (key, value) in D)
end


function Base.print(io::IO, prof::ObsProfile)
    TOML.print(io, struct_to_dict(prof))
end


function ObsProfile(filename)
    t = dict_to_tuple(TOML.parsefile(filename))
    return ObsProfile(;t...)
end






# Density methods


# Σ = h.values ./ (2π * log(10) * r .^ 2) # is equivalent because of derivative of log r
"""
    calc_Σ(log_r_bin, mass_per_annulus)

Calculate the surface density given the radii `log_r_bin` and the mass per annuli `mass_per_annulus`.
"""
function calc_Σ(log_r_bin, mass_per_annulus)
    r = 10 .^ log_r_bin
    As = π * diff(r .^ 2)

    Σ = mass_per_annulus ./ As 
	return Σ 
end


"""
    calc_Γ(log_rs, Σs)

Calculate the logarithmic slope of the density profile given the radii `log_rs` and the surface densities `Σs`.
"""
function calc_Γ(log_rs, Σs)
	d_log_r = lguys.gradient(log_rs)
	d_log_Σ = lguys.gradient(log10.(Σs))

	return d_log_Σ ./ d_log_r
end


function calc_M_in(log_r_bin, mass_per_annulus)
    M_in = cumsum(mass_per_annulus)
    return M_in 
end

function calc_Σ_mean(log_r_bin, M_in)
    r = 10 .^ log_r_bin[2:end]
	Areas = @. π * r^2
	Σ_bar = M_in ./ Areas
	return Σ_bar
end



function calc_Γ_max(Σ, Σ_m)
    Γ_max = @. 2*(1 - Σ / Σ_m)
    Γ_max[value.(Σ) .== 0] .= NaN
    return Γ_max
end


"""
    calc_properties(rs, r_units; weights=nothing, bins=20, normalization="mass")

Calculate the properties of a density profile given the radii `rs` and the units of the radii `r_units`.

"""
function calc_properties(rs; 
        r_units="", 
        weights=nothing, 
        bins=20, 
        normalization=true
    )

    if weights === nothing
        weights = ones(length(rs))
    end


    h = histogram(log10.(rs), bins, weights=weights, normalization=:count)
    log_r_bin = h.bins
    h.err[isnan.(h.err)] .= 0

    mass_per_annulus = h.values .± h.err
    counts = h.values

    if normalization
        mass_per_annulus = mass_per_annulus ./ sum(mass_per_annulus)
    end

    log_r = lguys.midpoint(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2
    log_r = log_r

    Σ = calc_Σ(log_r_bin, mass_per_annulus)

    M_in = calc_M_in(log_r_bin, mass_per_annulus)
    Σ_m = calc_Σ_mean(log_r_bin, M_in)

    Γ = calc_Γ(log_r, Σ)
    Γ_max = calc_Γ_max(Σ, Σ_m)


	log_Σ = log10.(Σ)
    log_Σ[Σ .== 0] .= NaN


    prof = ObsProfile(
        r_units = r_units,
        log_r = value.(log_r),
        log_r_bins = log_r_bin,
        counts = counts,
        M_in = value.(M_in),
        M_in_err = err.(M_in),
        mass_in_annulus = value.(mass_per_annulus),
        mass_in_annulus_err = err.(mass_per_annulus),
        Sigma = value.(Σ),
        Sigma_err = err.(Σ),
        Sigma_m = value.(Σ_m),
        Sigma_m_err = err.(Σ_m),
        log_Sigma = value.(log_Σ),
        log_Sigma_err = err.(log_Σ),
        Gamma = value.(Γ),
        Gamma_err = err.(Γ),
        Gamma_max = value.(Γ_max),
        Gamma_max_err = err.(Γ_max)
    )

    return prof
end



# RELLL


"""
    calc_r_ell(x, y, a, [b, ]PA)

computes the elliptical radius of a point (x, y) with respect to the center (0, 0) and the ellipse parameters (a, b, PA).
If using sky coordinates, x and y should be tangent coordinates.

Note that the position angle is the astronomy definition, i.e. measured from the North to the East (clockwise) in xi / eta.
"""
function calc_r_ell(x, y, args...)
    x_p, y_p = shear_points_to_ellipse(x, y, args...)

	r_sq = @. (x_p)^2 + (y_p)^2
	return sqrt.(r_sq)
end



"""
Transforms x and y into the sheared rotated frame of the ellipse.
"""
function shear_points_to_ellipse(x, y, a, b, PA)
	θ = @. deg2rad(90 - PA)
	x_p = @. x * cos(θ) + -y * sin(θ)
	y_p = @. x * sin(θ) + y * cos(θ)
    # scale
    x_p ./= a
    y_p ./= b

    return x_p, y_p
end

function shear_points_to_ellipse(x, y, ecc, PA)
    b = (1 - ecc^2)^(1/4)
    a = 1/b
    return shear_points_to_ellipse(x, y, a, b, PA)
end

"""
    calc_r_ell(x, y, ecc, PA)

Calculates the elliptical radius 
"""
function calc_r_ell(x, y, ecc, PA)
    b = (1 - ecc^2)^(1/4)
    a = 1/b
    println(ecc, " ", sqrt(1 - b^2 / a^2))
    return calc_r_ell(x, y, a, b, PA)
end



function calc_r_ell_sky(ra, dec, ecc, PA; kwargs...)
    b = (1 - ecc^2)^(1/4)
    a = 1/b
    return calc_r_ell_sky(ra, dec, a, b, PA; kwargs...)
end


"""
    calc_r_ell_sky(ra, dec, a, b, PA; weights=nothing, centre="mean", units="arcmin")

Given a set of sky coordinates (ra, dec), computes the elliptical radius of each point with respect to the centre of the ellipse defined by the parameters (a, b, PA).

Returns r_ell and the maximum radius within the convex hull of the points.
"""
function calc_r_ell_sky(ra, dec, a, b, PA; weights=nothing,
        centre="mean",
        units="arcmin"
    )
    if centre isa Tuple
        ra0, dec0 = centre
    else
        ra0, dec0 = calc_centre2D(ra, dec, centre, weights)
    end

    x, y = to_tangent(ra, dec, ra0, dec0)

    r_ell = calc_r_ell(x, y, a, b, PA)


    if units == "arcmin"
        r_ell = 60r_ell
    elseif units == "arcsec"
        r_ell = 3600r_ell
    elseif units == "deg"
        r_ell = 1r_ell
    else
        error("units not implemented: $units")
    end

    return r_ell
end


function calc_r_ell_sky(ra, dec; kwargs...)
    return calc_r_ell_sky(ra, dec, 0, 0; kwargs...)
end


function calc_r_max(ra, dec, args...; 
        centre="mean",
        units="arcmin",
        weights=nothing
    )

    if centre isa Tuple
        ra0, dec0 = centre
    else
        ra0, dec0 = calc_centre2D(ra, dec, centre, weights)
    end

    x, y = to_tangent(ra, dec, ra0, dec0)
    x_p, y_p = shear_points_to_ellipse(x, y, args...)
    
    hull = convex_hull(x_p, y_p)
    r_max = min_distance_to_polygon(hull...)
end



"""
    r_ell_max(xi, eta, ecc, PA)

returns the maximum radius of an ellipse within the convex boundary defined by the points.

Not implemented
"""
function r_ell_max(xi, eta, ecc, PA)
    error("not implemented")
end



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


"""
    distance_to_segment(a, b, p)

Distance from point `p` to the line segment defined by `a` and `b`.
all points are 2D vectors.
"""
function distance_to_segment(a, b, p=zeros(2))
    a = vec(a)
    b = vec(b)

    # work in origin at p
    a -= p
    b -= p

    # is the segment a point?
    l = norm(a - b)
    if l == 0
        return norm(a)  
    end

    # line unit vector
    n = (a - b) / l
    # projection along line
    t = dot(a, n) 

    if t < 0
        closest_point = a
    elseif t > l
        closest_point = b
    else
        closest_point = a - t * n
    end

    dist = norm(closest_point)

    return dist
end


function min_distance_to_polygon(x, y)
    min_dist = Inf
    N = length(x)
    for i in 1:N
        a = [x[i], y[i]]
        j = mod1(i + 1, N)
        b = [x[j], y[j]]

        dist = distance_to_segment(a, b)

        min_dist = min(min_dist, dist)
    end

    return min_dist
end

    



function mean_centre(ra, dec, mass)
	return mean(ra, weights(mass)), mean(dec, weights(mass))
end



function calc_centre2D(ra, dec, centre_method, weights=nothing)
    if weights === nothing
        weights = ones(length(ra))
    end
	if centre_method == "mean"
		ra0, dec0 = mean_centre(ra, dec, weights)
    else
        error("centre method not implemented: $centre_method")
	end

	return ra0, dec0
end

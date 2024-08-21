import Base: @kwdef

import LinearAlgebra: diag, dot, norm, normalize

using Measurements
import DensityEstimators: histogram
import StatsBase: weights, mean, StatsBase
import TOML


value(x::Measurement) = x.val
err(x::Measurement) = x.err
value(x) = x
err(x) = zero(x)


F = Float64


"""
An observed 2D density profile
"""
@kwdef mutable struct ObsProfile
    r_units::String

    log_r::Vector{F}
    log_r_bins::Vector{F}
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


function Base.print(io::IO, prof::ObsProfile)
    TOML.print(io, struct_to_dict(prof))
end

function ObsProfile(filename::String)
    t = dict_to_tuple(TOML.parsefile(filename))
    return ObsProfile(;t...)
end


"""
    ellipticity_to_aspect(ellipticity)

Converts the ellipticity to the aspect ratio (b/a) of an ellipse.
"""
function ellipticity_to_aspect(ellipticity)
    return 1 - ellipticity
end


"""
    calc_properties(rs, r_units; weights=nothing, bins=20, normalization="mass")

Calculate the properties of a density profile given the radii `rs` and the units of the radii `r_units`.


# Arguments
- `bins`: Passed to Arya.histogram. Bins in log r for histogram.
----------
rs : Vector{F}
    The radii of the profile.
r_units : String
    The units of the radii.
weights : Vector{F}, optional
    The weights of the radii.
bins : optional
    Bins passed to Arya.histogram.
normalization : String, optional
    The normalization of the profile.
    - :mass: normalizes the profile by the total mass.
    - :central: normalizes the profile by the mass within a central radius, r_centre
    - :none: no normalization.

"""
function calc_properties(rs; 
        r_units="", 
        weights=nothing, 
        bins=20, 
        normalization=:mass,
        r_centre=30,
    )

    if weights === nothing
        weights = ones(length(rs))
    end

    if any(rs .< 0)
        throw(DomainError("Radii must be positive."))
    end

    if length(rs) < 2
        throw(ArgumentError("Radii must have at least 2 elements."))
    end
    h = histogram(log10.(rs), bins, weights=weights, normalization=:none) # counting histograms
    log_r_bin = h.bins
    log_r = midpoints(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2

    h.err[isnan.(h.err)] .= 0

    mass_per_annulus = h.values .± h.err
    h_c = histogram(log10.(rs), bins, normalization=:none)
    counts = h_c.values

    if normalization == :mass
        mass_per_annulus = mass_per_annulus ./ sum(mass_per_annulus)
    elseif normalization == :central
        filt_cen = log_r .< log10(r_centre)
        Σ_m_cen = sum(mass_per_annulus[filt_cen]) ./ (π * r_centre^2)
        mass_per_annulus = mass_per_annulus ./ Σ_m_cen
    elseif normalization == :none
        mass_per_annulus = mass_per_annulus
    else
        error("normalization not implemented: $normalization")
    end

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


# Density methods
"""
    calc_Σ(log_r_bin, mass_per_annulus)

Calculate the surface density given the radii `log_r_bin` and the mass per annuli `mass_per_annulus`.
"""
function calc_Σ(log_r_bin, mass_per_annulus)
    r = 10. .^ log_r_bin
    As = π * diff(r .^ 2)

    # Σ = h.values ./ (2π * log(10) * r .^ 2) # is equivalent because of derivative of log r
    Σ = mass_per_annulus ./ As 
	return Σ 
end


"""
    calc_Γ(log_rs, Σs)

Calculate the logarithmic slope of the density profile given the radii `log_rs` and the surface densities `Σs`.
"""
function calc_Γ(log_rs, Σs)
	d_log_r = gradient(log_rs)
	d_log_Σ = gradient(log10.(Σs))

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



# Utility functions


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

function calc_r_ell_sky(ra, dec, ell, PA; kwargs...)
    aspect = ellipticity_to_aspect(ell)
    b = sqrt(aspect)
    a = 1/b
    return calc_r_ell_sky(ra, dec, a, b, PA; kwargs...)
end

function calc_r_ell_sky(ra, dec; kwargs...)
    return calc_r_ell_sky(ra, dec, 0, 0; kwargs...)
end


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
    shear_points_to_ellipse(x, y, a, b, PA)

Transforms x and y into the sheared rotated frame of the ellipse.
Position angle is measured from y axis in the direction of positive x.
"""
function shear_points_to_ellipse(x, y, a, b, PA)
    if a <= 0 || b <= 0
        throw(DomainError("a and b must be positive. Got a = $a, b = $b."))
    end

    θ = @. deg2rad(PA - 90)
	x_p = @. x * cos(θ) + -y * sin(θ)
	y_p = @. x * sin(θ) + y * cos(θ)
    # scale
    x_p = x_p ./ a
    y_p = y_p ./ b

    return x_p, y_p
end

function shear_points_to_ellipse(x, y, ell, PA)
    aspect = ellipticity_to_aspect(ell)
    b = sqrt(aspect)
    a = 1/b
    return shear_points_to_ellipse(x, y, a, b, PA)
end



"""
    calc_r_max(ra, dec, args...; centre="mean", weights=nothing)


"""
function calc_r_max(ra, dec, args...; 
        centre="mean",
        weights=nothing
    )

    ra0, dec0 = calc_centre2D(ra, dec, centre, weights)

    x, y = to_tangent(ra, dec, ra0, dec0)
    x_p, y_p = shear_points_to_ellipse(x, y, args...)

    if length(args) == 3
        a, b, _ = args
        aspect = b/a
        if aspect < 1
            aspect = 1/aspect
        end
    else
        aspect = ellipticity_to_aspect(args[1])
    end
    
    if isdefined(LilGuys, :convex_hull)
        hull = convex_hull(x_p, y_p)
        r_max = min_distance_to_polygon(hull...)
    else
        @warn "r_ell: Convex hull not defined. Using max radius. Load Polyhedra to enable convex hull."
        r_max = maximum(@. sqrt(x^2 + y^2)) ./ sqrt(aspect)
    end

    return r_max
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



function spherical_mean(ra, dec, weights=nothing)
    if length(ra) != length(dec)
        throw(DimensionMismatch("ra and dec must have the same length."))
    end

    pos = unit_vector(ra, dec)
    if weights === nothing
        mean_pos = centroid(pos)
    else
        mean_pos = centroid(pos, weights)
    end

    ra0, dec0, r= cartesian_to_sky(mean_pos...)
    if r == 0
        throw(DomainError("Mean position is at the origin."))
    end
    return ra0, dec0
end



"""
Calculates the 2d centre givin data.
"""
function calc_centre2D(ra, dec, centre_method, weights=nothing)
	if centre_method == "mean"
		ra0, dec0 = spherical_mean(ra, dec, weights)
    elseif centre_method isa Tuple
        ra0, dec0 = centre_method
    else
        error("centre method not implemented: $centre_method")
	end

	return ra0, dec0
end



"""
    to_orbit_coords(ra, dec, ra0, dec0, PA) 

Given the position angle of an orbit vector at (ra0, dec0) calculates a rotated
sky frame centred on RA, DEC and with the x-axis aligned with the orbit.
The orbital position angle can be found in the orbit properties file from
analyze_orbit.jl notebook.
"""
function to_orbit_coords(ra, dec, ra0::Real, dec0::Real, PA::Real)
	# want to rotate to dec, ra of centre, then rotate 
	α = deg2rad(ra0)
	δ = deg2rad(dec0)
	ϖ = deg2rad(90 - PA)
	Rmat = Rx_mat(ϖ) * Ry_mat(δ) * Rz_mat(-α) 

	coords = unit_vector(ra, dec)
	coords =  Rmat * coords
	ra, dec, _ = cartesian_to_sky(coords[1, :], coords[2, :], coords[3, :])

	ra .-= 360 * (ra .> 180)

	ra, dec
end

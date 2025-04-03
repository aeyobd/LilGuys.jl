import Base: @kwdef, filter

import LinearAlgebra: diag, dot, norm, normalize, ⋅, × 

import TOML



F = Float64


"""
An observed 2D density profile

"""
@kwdef mutable struct StellarDensityProfile
    R_units::String
    mass_units::String = ""

    log_R::Vector{F}
    log_R_bins::Vector{F}
    counts::Vector{F} = F[]

    log_Sigma::Vector{Measurement{F}}

    Gamma::Vector{Measurement{F}} = Measurement{F}[]

    log_m_scale::F = 0.
    "log radius shift used for normalization"
    log_R_scale::F = 0.

    annotations::Dict{String, Any} = Dict{String, Any}()
end


"""
An observed 2D mass profile (cumulative)
"""
@kwdef mutable struct StellarMassProfile
    R_units::String
    mass_units::String = ""

    log_R::Vector{F}
    M_in::Vector{F}
    M_in_err::Vector{F}

    Sigma_m::Vector{F} = []
    Sigma_m_err::Vector{F} = []

    log_m_scale::F = 0.
    "log radius shift used for normalization"
    log_R_scale::F = 0.

    annotations::Dict = Dict()
end



"""
    StellarDensityProfile(radii; arguments...)

Calculate the properties of a density profile given the radii `Rs` and the units of the radii `R_units`.


# Arguments
- `weights`: The weights of the radii.
- `bins`: Passed to Arya.histogram. Bins in log r for histogram.
- `errors`: Passed to Arya.histogram. Histogram error method. May be :bernoulli, :poisson, :weighted, or :none. Defaults to :weighted.
- `normalization`: The normalization of the profile.
    - :mass: normalizes the profile by the total mass.
    - :central: normalizes the profile by the mass within a central radius, R_centre
    - :none: no normalization.
- `bins_centre`: The number of bins to use for centre normalization.
- `distance`: distance (optional)
- `sigma_v`: velocity dispersion (just as an annotation)
- `R_units`: The units of the radii, entirely for self-documentation currently.
"""
function StellarDensityProfile(Rs; 
        weights=nothing, 
        bins=nothing, 
        normalization=:none,
        bins_centre=3,
        R_units="", 
        distance=NaN,
        sigma_v=NaN,
        quantiles = [0.01, 0.02, 0.05, 0.10, 0.16, 0.25, 0.5, 0.75, 0.84, 0.9, 0.95, 0.98, 0.99],
        errors=:weighted,
        kwargs...
    )

    if weights === nothing
        weights = ones(length(Rs))
    end

    if any(Rs .< 0)
        throw(DomainError("Radii must be positive."))
    end

    if length(Rs) < 2
        throw(ArgumentError("Radii must have at least 2 elements."))
    end

    Rs, weights = Interface.filter_nans(Rs, weights)

    log_R_bin, mass_per_annulus, mass_per_annulus_err = histogram(
        log10.(Rs), bins, 
        weights=weights, normalization=:none, errors=errors
    )
    _, counts, _ = histogram(log10.(Rs), log_R_bin, normalization=:none)


    # set nans one for now, revert later
    mass_per_annulus_err[counts .== 0] .= 1

    log_R_m = midpoints(log_R_bin)
    δ_log_R = -diff(log_R_bin) / 2

    Σ, Σ_err = calc_Σ_from_hist(log_R_bin, mass_per_annulus, mass_per_annulus_err)

    log_Σ = log10.(Σ)
    log_Σ[Σ .== 0] .= -Inf
    log_Σ_em = log_Σ .- log10.(max.(Σ .- Σ_err, 0))
    log_Σ_ep = log10.(Σ .+ Σ_err) .- log_Σ

    log_Sigma = Measurement.(log_Σ, log_Σ_em, log_Σ_ep)

    Γ = calc_Γ(Measurement.(log_R_m, δ_log_R), log_Sigma)

    prof = StellarDensityProfile(;
        R_units = R_units,
        log_R = log_R_m,
        log_R_bins = log_R_bin,
        counts = counts,
        log_Sigma = log_Sigma,
        Gamma = Γ,
        annotations = Dict(
            "normalization" => string(normalization),
           ),
    )

    prof = normalize(prof, mass_per_annulus, normalization, bins_centre=bins_centre)
    return prof
end


"""
    StellarDensityProfile(snap::Snapshot; R_units="kpc", x_vec=[1, 0, 0], y_vec=[0, 1, 0], kwargs...)

Creates a StellarDensityProfile from a snapshot. 
If the units is set to kpc, than the profile is calculated in the xy plane defined
by the vectors x_vec and y_vec. Otherwise, the snapshot is projected as is to the sky by `to_gaia` and then the profile is calculated using the circular radius
of the on-sky tangent plane. Kwarguments are passed to the StellarDensityProfile constructor for a list of radii.
"""
function StellarDensityProfile(snap::Snapshot;
        R_units = "kpc",
        x_vec = [1, 0, 0],
        y_vec = [0, 1, 0],
        kwargs...
    )

    if R_units == "kpc"
        N = size(snap.positions, 2)
        ra = [dot(snap.positions[:, i], x_vec) for i in 1:N]
        dec = [dot(snap.positions[:, i], y_vec) for i in 1:N]

        weights = snap.weights
        ra0 = snap.x_cen ⋅ x_vec
        dec0 = snap.x_cen ⋅ y_vec
    elseif R_units ∈ ["arcmin", "arcsec", "deg"]
        projected = to_gaia(snap, add_centre=true) # centre has zero weight
        ra = projected.ra
        dec = projected.dec
        ra0 = projected.ra[1]
        dec0 = projected.dec[1]
        weights = projected.weights
    else
        throw(ArgumentError("R_units not implemented: $R_units"))
    end

    xi, eta = to_tangent(ra, dec, ra0, dec0)
    r = sqrt.(xi .^ 2 + eta .^ 2)

    annotations = Dict("time" => snap.time)
    if :annotations ∈ [k for (k, v) in kwargs]
        for (key, val) in kwargs[:annotations]
            annotations[string(key)] = val
        end
    end

    return StellarDensityProfile(r; R_units=R_units, weights=weights, annotations=annotations)
end



function Base.print(io::IO, prof::StellarDensityProfile)
    TOML.print(io, struct_to_dict(prof))
end


"""
    StellarDensityProfile(filename::String; kwargs...)

Reads a StellarDensityProfile from a TOML file.
"""
function StellarDensityProfile(filename::String; kwargs...)
    t = TOML.parsefile(filename)
    t = merge(t, kwargs)
    t = collapse_errors(t)
    t = dict_to_tuple(t)

    return StellarDensityProfile(;t...)
end





"""
    nan_uncertainty(v)

Returns the uncertainty of a value `v` with zero uncertainties as NaN
"""
function nan_uncertainty(v::AbstractVector{<:Real})
    return fill(NaN, length(v))
end


"""
    ellipticity_to_aspect(ellipticity)

Converts the ellipticity to the aspect ratio (b/a) of an ellipse.
"""
function ellipticity_to_aspect(ellipticity::Real)
    return 1 - ellipticity
end



# Density methods
"""
    calc_Σ_from_hist(log_R_bin, mass_per_annulus)

Calculate the surface density given the radii `log_R_bin` and the mass per annuli `mass_per_annulus`.
"""
function calc_Σ_from_hist(log_R_bin::AbstractVector{<:Real}, mass_per_annulus::AbstractVector{<:Real}, mass_per_annulus_err)
    r = 10. .^ log_R_bin
    As = π * diff(r .^ 2)

    Σ = mass_per_annulus ./ As 
    Σ_err = mass_per_annulus_err ./ As
    return Σ, Σ_err
end


"""
    calc_Σ_from_1D_density(log_R_bin, density1d)

Calculate the surface density given the radii `log_R_bin` and the 1D density `density1d`.
"""
function calc_Σ_from_1D_density(log_R_bin::AbstractVector{<:Real}, density1d::AbstractVector{<:Real})
    r = 10 .^ log_R_bin
    As = π * diff(r .^ 2)

    Σ = density ./ (2π * log(10) * r .^ 2)
    return Σ
end


"""
    calc_Γ(log_Rs, Σs)

Calculate the logarithmic slope of the density profile given the radii `log_Rs` and the surface densities `Σs`.
"""
function calc_Γ(log_Rs::AbstractVector{<:Real}, log_Σ::AbstractVector{<:Measurement})
    d_log_R = gradient(log_Rs)
    d_log_Σ = gradient(log_Σ)

    return d_log_Σ ./ d_log_R
end


"""
    calc_M_partial(log_Rs, log_R_bins, log_R_min; weights)

Computes the amount of mass between each bin and the next bin midpoints
"""
function calc_M_partial(log_Rs, log_R_bins, log_R_m; weights=ones(length(log_Rs)))
    M_partials = zeros(length(log_R_m))
    M_errs = zeros(length(log_R_m))

    for i in eachindex(M_partials)
        filt = log_Rs .< log_R_m[i]
        filt .&= log_Rs .>= log_R_bins[i]
        M_partials[i] = sum(weights[filt])
        M_errs[i] = sqrt.(sum(weights[filt] .^ 2))
    end

    return M_partials .± M_errs
end


"""
    calc_Σ_mean_from_hist(log_R, M_in)

Calculates the mean surface density from the interior mass to a radius
"""
function calc_Σ_mean_from_hist(log_R::AbstractVector{<:Real}, M_in::AbstractVector{<:Real})
    R = 10 .^ log_R
    Areas = @. π*R^2
    Σ_bar = M_in ./ Areas
    return Σ_bar
end


"""
Given the surface density Σ and the mean surface density Σ_m, calculates the maximum slope of the density profile.
""" 
function calc_Γ_max(Σ, Σ_m)
    Γ_max = @. 2*(1 - Σ / Σ_m)
    Γ_max[Σ .== 0] .= NaN
    return Γ_max
end



# Utility functions


"""
    calc_R_ell_sky(ra, dec, a, b, PA; weights=nothing, centre="mean", units="arcmin")

Given a set of sky coordinates (ra, dec), computes the elliptical radius of each point with respect to the centre of the ellipse defined by the parameters (a, b, PA).
"""
function calc_R_ell_sky(ra, dec, a, b, PA; weights=nothing,
        centre="mean",
        units="arcmin"
    )
    ra0, dec0 = calc_centre2D(ra, dec, centre, weights)

    x, y = to_tangent(ra, dec, ra0, dec0)

    R_ell = calc_R_ell(x, y, a, b, PA)


    if units == "arcmin"
        R_ell = 60R_ell
    elseif units == "arcsec"
        R_ell = 3600R_ell
    elseif units == "deg"
        R_ell = 1R_ell
    else
        error("units not implemented: $units")
    end

    return R_ell
end


function calc_R_ell_sky(ra, dec, ell, PA; kwargs...)
    aspect = ellipticity_to_aspect(ell)
    b = sqrt(aspect)
    a = 1/b
    return calc_R_ell_sky(ra, dec, a, b, PA; kwargs...)
end


function calc_R_ell_sky(ra, dec; kwargs...)
    return calc_R_ell_sky(ra, dec, 0, 0; kwargs...)
end


"""
    calc_R_ell(x, y, a, [b, ]PA)

computes the elliptical radius of a point (x, y) with respect to the center (0, 0) and the ellipse parameters (a, b, PA).
If using sky coordinates, x and y should be tangent coordinates.

Note that the position angle is the astronomy definition, i.e. measured from the North to the East (clockwise) in xi / eta.
"""
function calc_R_ell(x, y, args...)
    x_p, y_p = shear_points_to_ellipse(x, y, args...)

    R_sq = @. (x_p)^2 + (y_p)^2
    return sqrt.(R_sq)
end



"""
    shear_points_to_ellipse(x, y, a, b, PA)

Transforms x and y into the sheared rotated frame of the ellipse.
Position angle is measured from y axis in the direction of positive x.
(North to East on sky).
"""
function shear_points_to_ellipse(x, y, a, b, PA)
    @assert_same_size x y

    if a <= 0 || b <= 0
        throw(DomainError("a and b must be positive. Got a = $a, b = $b."))
    end

    # rotate
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
    spherical_mean(ra, dec, weights=nothing)

Calculates the spherical mean of a set of sky coordinates (ra, dec) with optional weights. The mean is calculated by taking the average of the unit vectors of the points, so avoids problems with periodicity and spherical geometry.
"""
function spherical_mean(ra::AbstractVector{<:Real}, dec::AbstractVector{<:Real}, weights::Union{Nothing, AbstractVector{<:Real}}=nothing)
    if length(ra) != length(dec)
        throw(DimensionMismatch("ra and dec must have the same length."))
    end

    pos = unit_vector(ra, dec)
    if weights === nothing
        mean_pos = centroid(pos)
    else
        mean_pos = centroid(pos, weights)
    end

    ra0, dec0, r = cartesian_to_sky(mean_pos...)
    if r == 0
        throw(DomainError("Mean position is at the origin."))
    end
    return ra0, dec0
end



"""
    calc_centre2D(ra, dec, centre_method, weights=nothing)

Calculates the centre of vectors of ra and dec (with weights) using one of the
following methods:
- "mean": calculates the mean of the vectors
- Tuple{ra0, dec0}: uses the given ra0 and dec0 as the centre (bypasses calculation)
This is a utility function.

"""
function calc_centre2D(ra::AbstractVector{<:Real}, dec::AbstractVector{<:Real}, centre_method, weights=nothing)
    if centre_method == "mean"
        ra0, dec0 = spherical_mean(ra, dec, weights)
    elseif centre_method isa Tuple
        ra0, dec0 = centre_method
    else
        throw(ArgumentError("centre method not implemented: $centre_method"))
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
	Rmat = Rx_mat(-ϖ) * Ry_mat(δ) * Rz_mat(-α) 

	coords = unit_vector(ra, dec)
	coords =  Rmat * coords
	ra, dec, _ = cartesian_to_sky(coords[1, :], coords[2, :], coords[3, :])

	ra .-= 360 * (ra .> 180)

	ra, dec
end



function to_orbit_coords(ra::Real, dec::Real, ra0::Real, dec0::Real, PA::Real)
    ra, dec =  to_orbit_coords([ra], [dec], ra0, dec0, PA)
    return ra[1], dec[1]
end


"""
    normalize(prof::StellarDensityProfile, normalization=:none)

Returns a normalized version of the stellar profile.
Options are
- `:mass` normalizes profile by total mass

"""
function normalize(prof::StellarDensityProfile, mass_per_annulus, normalization=:none; bins_centre=3)
    if normalization isa String
        normalization = Symbol(normalization)
    end

    @assert prof.log_m_scale == 0. "can only normalize profile once"

    if normalization == :mass
        m_scale = 1/ sum(mass_per_annulus)
    elseif normalization == :central
        filt_cen = 1:bins_centre
        R_centre = 10 .^ prof.log_R_bins[bins_centre+1]

        Σ_m_cen = sum(mass_per_annulus[filt_cen]) ./ (π * R_centre^2)
        m_scale = 1 / Σ_m_cen
    elseif normalization == :none
        m_scale = 1
    else
        error("normalization not implemented: $normalization")
    end

    return scale(prof, 1.0, m_scale; _normalization=string(normalization))
end

"""
    scale(prof::StellarDensityProfile, R_scale::Real, m_scale::Real)

Scales the profile by a factor of `R_scale` in radius and `m_scale` in mass,
returning a new profile.
"""
function scale(prof::StellarDensityProfile, R_scale::Real, m_scale::Real; _normalization=nothing)
    prof_new = deepcopy(prof)
    prof_new.log_m_scale += log10(m_scale)
    prof_new.log_R_scale += log10(R_scale)
    prof_new.log_R .+= log10(R_scale)
    prof_new.log_R_bins .+= log10(R_scale)
    prof_new.log_Sigma .+= log10(m_scale) .- 2log10(R_scale)
    if !isnothing(_normalization)
        prof_new.annotations["normalization"] = _normalization
    end

    return prof_new
end

"""
filters a profile by the given bitarray (indexed by bin)
"""
function filter_by_bin(prof::StellarDensityProfile, bin_selection)
    edge_filt = edge_from_midpoint_filter(bin_selection)
    filt = bin_selection

    prof_new = deepcopy(prof)
    prof_new.log_R = prof.log_R[filt]
    prof_new.log_R_bins = prof.log_R_bins[edge_filt]
    prof_new.counts = prof.counts[filt]
    prof_new.log_Sigma = prof.log_Sigma[filt]

    if length(prof.Gamma) > 0
        prof_new.Gamma = prof.Gamma[filt]
    end

    return prof_new
end



function filter_empty_bins(prof::StellarDensityProfile)
    idxs = find_longest_consecutive_finite(prof.log_Sigma)

    filt = eachindex(prof.log_R) .∈ Ref(idxs)
    return filter_by_bin(prof, filt)
end


"returns the bin edges if x is a bitarray filter on bins"
function edge_from_midpoint_filter(x)
    return [false; x] .| [x; false]
end


"""
    find_longest_consequtive_finite(x)

Returns the longest consecutive sequence of finite elements in x
as a int-range.
"""
function find_longest_consecutive_finite(x)
    max_len = 0
    max_start = 0
    max_end = 0
    current_len = 0
    current_start = 0

    for i in eachindex(x)
        if isfinite(x[i])
            if current_len == 0
                current_start = i
            end
            current_len += 1
        else
            if current_len > max_len
                max_len = current_len
                max_start = current_start
                max_end = i - 1
            end
            current_len = 0
        end
    end

    # Check the last segment after loop ends
    if current_len > max_len
        max_len = current_len
        max_start = current_start
        max_end = lastindex(x)
    end

    return max_len > 0 ? (max_start:max_end) : nothing
end

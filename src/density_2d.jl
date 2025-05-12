import Base: @kwdef, filter
import TOML
import LinearAlgebra: diag, dot, norm, normalize, ⋅, × 


"""
An observed 2D stellar density profile

"""
@kwdef mutable struct SurfaceDensityProfile
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


units(prof::SurfaceDensityProfile) = (length => prof.R_units, mass => prof.mass_units)
log_radii(prof::SurfaceDensityProfile) = prof.log_R
radii(prof::SurfaceDensityProfile) = 10 .^ prof.log_R
log_surface_density(prof::SurfaceDensityProfile) = prof.log_Sigma
surface_density(prof::SurfaceDensityProfile) = 10 .^prof.log_Sigma
log_surface_density_err(prof::SurfaceDensityProfile) = sym_error.(prof.log_Sigma)
surface_density_err(prof::SurfaceDensityProfile) = sym_error.(density_2D(prof))



"""
An observed 2D mass profile (cumulative)
"""
@kwdef mutable struct CylMassProfile
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
    SurfaceDensityProfile(radii; arguments...)

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
- `R_units`: The units of the radii, entirely for self-documentation currently.
- `annotations`: Additional notes to be added to the profile.
"""
function SurfaceDensityProfile(Rs; 
        weights=nothing, 
        bins=nothing, 
        normalization=:none,
        bins_centre=3,
        R_units="", 
        distance=NaN,
        errors=:weighted,
        annotations=Dict{String, Any}(),
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

    annotations["normalization"] = string(normalization)

    prof = SurfaceDensityProfile(;
        R_units = R_units,
        log_R = log_R_m,
        log_R_bins = log_R_bin,
        counts = counts,
        log_Sigma = log_Sigma,
        Gamma = Γ,
        annotations = annotations
    )

    prof = normalize(prof, mass_per_annulus, normalization, bins_centre=bins_centre)
    return prof
end


"""
    SurfaceDensityProfile(snap::Snapshot; R_units="kpc", x_vec=[1, 0, 0], y_vec=[0, 1, 0], kwargs...)

Creates a SurfaceDensityProfile from a snapshot. 
If the units is set to kpc, than the profile is calculated in the xy plane defined
by the vectors x_vec and y_vec. Otherwise, the snapshot is projected as is to the sky by `to_gaia` and then the profile is calculated using the circular radius
of the on-sky tangent plane. Kwarguments are passed to the SurfaceDensityProfile constructor for a list of radii.
"""
function SurfaceDensityProfile(snap::Snapshot;
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

    annotations = Dict{String, Any}("time" => snap.time)
    if :annotations ∈ keys(kwargs)
        for (key, val) in kwargs[:annotations]
            annotations[string(key)] = val
        end
    end

    return SurfaceDensityProfile(r; R_units=R_units, weights=weights, annotations=annotations)
end



function Base.print(io::IO, prof::SurfaceDensityProfile)
    TOML.print(io, struct_to_dict(prof))
end


"""
    SurfaceDensityProfile(filename::String; kwargs...)

Reads a SurfaceDensityProfile from a TOML file.
"""
function SurfaceDensityProfile(filename::String; kwargs...)
    t = TOML.parsefile(filename)
    t = merge(t, Dict(string(k) => v for (k, v) in kwargs))
    t = Dict{String, Any}(t)
    t = collapse_errors(t)
    t = dict_to_tuple(t)

    return SurfaceDensityProfile(;t...)
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


# not presently used
# """
#     calc_Σ_from_1D_density(log_R_bin, density1d)
# 
# Calculate the surface density given the radii `log_R_bin` and the 1D density `density1d`.
# """
# function calc_Σ_from_1D_density(log_R_bin::AbstractVector{<:Real}, density1d::AbstractVector{<:Real})
#     r = 10 .^ log_R_bin
# 
#     Σ = density ./ (2π * log(10) * r .^ 2)
#     return Σ
# end
# 

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



"""
    normalize(prof::SurfaceDensityProfile, normalization=:none)

Returns a normalized version of the stellar profile.
Options are
- `:mass` normalizes profile by total mass

"""
function normalize(prof::SurfaceDensityProfile, mass_per_annulus, normalization=:none; bins_centre=3)
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
    scale(prof::SurfaceDensityProfile, R_scale::Real, m_scale::Real)

Scales the profile by a factor of `R_scale` in radius and `m_scale` in mass,
returning a new profile.
"""
function scale(prof::SurfaceDensityProfile, R_scale::Real, m_scale::Real; _normalization=nothing)
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
function filter_by_bin(prof::SurfaceDensityProfile, bin_selection::UnitRange)
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



function filter_empty_bins(prof::SurfaceDensityProfile)
    idxs = find_longest_consecutive_finite(log_density_2D(prof))

    return filter_by_bin(prof, idxs)
end


"returns the bin edges if x is a bitarray filter on bins"
function edge_from_midpoint_filter(x::UnitRange)
    return (x.start):(x.stop+1)
end


"""
    find_longest_consequtive_finite(x)

Returns the longest consecutive sequence of finite elements in x
as a int-range.
Returns nothing if there are no finite elements in x.
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
    
    if max_len > 0
        return max_start:max_end
    else
        return 1:0 # null selection
    end
end

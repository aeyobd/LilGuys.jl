import TOML
import LinearAlgebra: ×, ⋅


"""
A struct representing a 3-dimensional 
mass profile (likely of a simulation snapshot).
All properties are in code units.
"""
@kwdef struct DensityProfile
    "Log radius of bin middles"
    log_r::Vector{F}

    "Log radius of bin edges"
    log_r_bins::Vector{F}

    "Number of particles in each bin"
    counts::Vector{Int}

    "Effective sample size in bin (if weighted)"
    ess::Vector{F} = F[]

    "Density in each bin"
    rho::Vector{Measurement{F}}

    "3D velocity dispersion"
    sigma_v::Vector{F} = []

    "anisotropy parameter"
    beta::Vector{F} = []

    "snapshot time"
    time::F = NaN

    "Additional annotations"
    annotations::Dict{String, Any} = Dict{String, Any}()
end

"the radii representing the middle of each bin for the profile"
log_radii(prof::DensityProfile) = prof.log_r
radii(prof::DensityProfile) = 10 .^ log_radii(prof)

log_radius_bins(prof::DensityProfile) = prof.log_r_bins
radius_bins(prof::DensityProfile) = 10 .^ log_radius_bins(prof)

counts_per_bin(prof::DensityProfile) = prof.counts

densities(prof::DensityProfile) = prof.rho
densities_err(prof::DensityProfile) = maximum.(error_interval.(prof.rho))
log_densities(prof::DensityProfile) = log10.(densities(prof))
log_densities_err(prof::DensityProfile) = maximum.(error_interval.(log10.(densities(prof))))

velocity_dispersion_3D(prof::DensityProfile) = prof.sigma_v
anisotropies(prof::DensityProfile) = prof.beta
effective_sample_sizes(prof::DensityProfile) = prof.ess


function Base.print(io::IO, prof::DensityProfile)
    TOML.print(io, struct_to_dict(prof))
end


function DensityProfile(filename::String)
    t = dict_to_tuple(TOML.parsefile(filename))
    return DensityProfile(;t...)
end


"""
    DensityProfile(snap; bins=100, filt_bound=:simple)

Given a snapshot, computes the density, circular velocity, and energy; 
returning an DensityProfileobject.

`bins` may be an integer, array, or function and is passed to `histogram`

Note: does not include binning error currently.
"""
function DensityProfile(snap::Snapshot;
        bins=nothing,
        filt_bound=:simple,
    )

    if filt_bound != :false
        filt = bound_particles(snap, method=filt_bound)
        snap = snap[filt]
    end

    if length(snap) == 0
        throw(ArgumentError("No bound particles in snapshot"))
    end

    log_r_snap = log10.(radii(snap))

    bins, hist, err = histogram(log_r_snap, bins, weights=snap.masses, errors=:weighted)
    log_r_bins = bins
    mass_in_shell = hist 
    mass_in_shell_err  = err

    _, counts, _ = histogram(log_r_snap, bins)

    log_r = midpoints(log_r_bins)
    r_bins = 10 .^ log_r_bins
    rel_err = mass_in_shell_err ./ mass_in_shell

    rho = density_from_hist(r_bins, mass_in_shell)
    rho_err = rho .* rel_err 

    return DensityProfile(
        log_r=log_r,
        log_r_bins=log_r_bins,
        counts=counts,
        rho=Measurement{F}.(rho, rho_err),
    )
end


"""
    DensityProfile(snap, weights; bins, filt_bound)

Weighted 3D density profile.
"""
function DensityProfile(snap::Snapshot, weights::AbstractVector{<:Real};
        bins=nothing,
        filt_bound=false,
    )

    if filt_bound
        filt = bound_particles(snap)
        snap = snap[filt]
    end

    if length(snap) == 0
        throw(ArgumentError("No bound particles in snapshot"))
    end

    log_r_snap = log10.(radii(snap))

    bins, hist, err = histogram(log_r_snap, bins, weights=weights, errors=:weighted)
    log_r_bins = bins
    mass_in_shell = hist 
    mass_in_shell_err  = err

    Nb = length(mass_in_shell)
    ess = zeros(Nb)
    counts = zeros(Nb)
    for i in 1:length(Nb)
        filt = log_r_snap .>= log_r_bins[i]
        filt .&= log_r_snap .< log_r_bins[i+1]
        w = weights[filt]
        counts[i] = sum(filt)
        ess[i] = effective_sample_size(w)
    end

    log_r = midpoints(log_r_bins)
    r_bins = 10 .^ log_r_bins
    rel_err = mass_in_shell_err ./ mass_in_shell

    rho = density_from_hist(r_bins, mass_in_shell)
    rho_err = rho .* rel_err 

    return DensityProfile(
        log_r=log_r,
        log_r_bins=log_r_bins,
        counts=counts,
        rho=Measurement{F}.(rho, rho_err),
        ess = ess,
        time = snap.time
    )
end


"""
    density_from_hist(r_bins, counts)
"""
function density_from_hist(r_bins::AbstractVector{<:Real}, counts::AbstractVector{<:Real}) 
    if length(r_bins) != length(counts) + 1
        throw(DimensionMismatch("Bins must be one longer than counts, got sizes $(length(r_bins)), $(length(counts))"))
    end

    volumes = 4π/3 * diff(r_bins .^ 3)
    return counts ./ volumes
end



""" 
	radial_velocities(snap)

returns the radial velocities relative to the snapshot centre in code units
"""
function radial_velocities(snap)
    return radial_velocities(snap.positions, snap.velocities, x_cen=snap.x_cen, v_cen=snap.v_cen)
end


"""
    radial_velocities(positions, velocities; x_cen=zeros(3), v_cen=zeros(3))

Calculates the radial velocities relative to x_cen, v_cen.
"""
function radial_velocities(positions::AbstractMatrix{<:Real}, velocities::AbstractMatrix{<:Real}; x_cen=zeros(3), v_cen=zeros(3))
    @assert_same_size positions velocities
    @assert_3vector positions

    x_vec = positions .- x_cen
    v_vec = velocities .- v_cen

    # normalize
    x_hat = x_vec ./ radii(x_vec)'

    # dot product
    v_r = sum(x_hat .* v_vec, dims=1)

    # matrix -> vector
    v_r = dropdims(v_r, dims=1)
    
    return v_r 
end



"""
    β_prof(snap; r_bins)

Compute the velocity dispersion and anisotropy profiles
for the snapshot assuming spherical symmetry.
Return a tuple of the 3D velocity dispersion and anisotropy in each bin
"""
function β_prof(snap; r_bins)
    vel = snap.velocities .- snap.v_cen
    pos = snap.positions .- snap.x_cen

    rs = radii(snap)
    v_r, v_θ, v_ϕ = to_spherical_velocities(pos, vel, rs)

    Nb = length(r_bins) - 1
    βs = zeros(Nb)
    σ3d = zeros(Nb)

    for i in 1:Nb
        filt = rs .> r_bins[i]
        filt .&= rs .<= r_bins[i+1]

        σ2_r = std(v_r[filt])^2
        σ2_t = std(v_θ[filt])^2 + std(v_ϕ[filt])^2

        βs[i] = 1 - σ2_t / 2σ2_r
        σ3d[i] = σ2_t + σ2_r
    end

    return σ3d, βs
end



"""
    spherical_unit_vectors(pos::AbstractVector{<:Real])

Compute the spherical unit vectors at the given position. Return 
`r_hat`, `θ_hat`, and `ϕ_hat` as a tuple
"""
function spherical_unit_vectors(pos::AbstractVector{<:Real}, rs=radii(pos))
    r_hat = pos ./ rs

    ϕ_vec = r_hat × [0,0,1]
    ϕ_norm = radii(ϕ_vec)
    if ϕ_norm ≈ 0
        ϕ_hat = zeros(3)
    else
        ϕ_hat = ϕ_vec ./ ϕ_norm
    end

    θ_vec = r_hat × ϕ_hat # should be perpendicular so 
    @assert radii(θ_vec) ≈ 1.0
    θ_hat = θ_vec ./ radii(θ_vec)

    return r_hat, θ_hat, ϕ_hat
end


"""
    to_spherical_velocities(positions, velocities[, rs])

Return the spherical velocities (r, θ, ϕ) given a matrix of positions
and velocities. The norm of eash position may be precomputed.
"""
function to_spherical_velocities(positions, velocities, rs=radii(positions))
    N = size(positions, 2)

    v_r = zeros(N)
    v_θ = zeros(N)
    v_ϕ = zeros(N)

    for i in 1:N
        pos = positions[:, i]
        vel = velocities[:, i]
        r_hat, θ_hat, ϕ_hat = spherical_unit_vectors(pos, rs[i])

        v_r[i] = vel ⋅ r_hat
        v_θ[i] = vel ⋅ θ_hat
        v_ϕ[i] = vel ⋅ ϕ_hat
    end

    return v_r, v_θ, v_ϕ
end



"""
    scale(prof::DensityProfile3D, r_scale, m_scale, m_scale_pot)

Scale the profile by the given factors, returning a new
"""
function scale(prof::DensityProfile, r_scale::Real,  m_scale::Real, m_scale_pot::Real)
    ρ_scale = m_scale / r_scale^3
    v_scale = sqrt(m_scale_pot / r_scale)


    return DensityProfile(
        log_r = prof.log_r .+ log10(r_scale),
        log_r_bins = prof.log_r_bins .+ log10(r_scale),
        rho = prof.rho .* ρ_scale,
        counts = prof.counts,
        ess = prof.ess,
        sigma_v = prof.sigma_v .* v_scale,
        beta = prof.beta,
        time = prof.time,
        annotations = prof.annotations
   )
end

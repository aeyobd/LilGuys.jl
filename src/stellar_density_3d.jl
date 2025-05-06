"""
    StellarDensity3D

Represents a stellar density profile in 3D,
"""
@kwdef struct StellarDensity3D
    "the log of the radius of the bin"
    log_r::Vector{F}

    "the log radii used for bins"
    log_r_bins::Vector{F}

    "the stellar density"
    rho::Vector{Measurement{F}}

    "number of particles in bin"
    counts::Vector{F}

    "effective sample size of bin (accounting for weights)"
    ess::Vector{F}

    "snapshot timestep"
    time::F = NaN

    annotations::Dict{String, Any} = Dict{String, Any}()
end


"""
    StellarDensity3D(snap; log_r_bins)

Calculate the 3D stellar density profile of the snapshot assuming the provided
bins. Bins passed to `histogram` and errors assumed to be `:weighted`
"""
function StellarDensity3D(snap;  bins=nothing)

    log_r_snap = log10.(radii(snap))

    bins, hist, err = histogram(log_r_snap, bins, weights=snap.weights, errors=:weighted)
    log_r_bins = bins
    mass_in_shell = hist 
    mass_in_shell_err  = err

    Nb = length(mass_in_shell)
    counts = zeros(Nb)
    ess = zeros(Nb)

    for i in 1:length(Nb)
        filt = log_r_snap .>= log_r_bins[i]
        filt .&= log_r_snap .< log_r_bins[i+1]
        w = snap.weights[filt]
        counts[i] = sum(filt)
        ess[i] = effective_sample_size(w)
    end
        

    log_r = midpoints(log_r_bins)
    r_bins = 10 .^ log_r_bins
    r = 10 .^ log_r
    rel_err = mass_in_shell_err ./ mass_in_shell

    rho = density_from_hist(r_bins, mass_in_shell)
    rho_err = rho .* rel_err 

    return StellarDensity3D(
        log_r=log_r,
        log_r_bins=log_r_bins,
        rho = Measurement{F}.(rho, rho_err),
        counts = counts,
        ess = ess,
        time = snap.time,
   )
end



"A structure containing basic (scalar) properties of a snapshot"
@kwdef struct StellarScalars
    "time since last peri"
    delta_t::F

    "break radius alla peñarrubia"
    r_break::F

    "1d velocity dispersion"
    sigma_v::F

    "3d radius cutoff for sigma_v"
    r_max_sigma::F

    "total bound stellar mass"
    bound_mass::F

    "quantiles used for `r_quantile`"
    quantiles::Vector{F}

    "radius containing `quantile` of stellar mass"
    r_quantiles::Vector{F}

    "snapshot time"
    time::F
end


"""
    StellarScalars(snap::Snapshot)

Compute essential scalar properties of the snapshot.

- `delta_t::Real=NaN`: The time since the last pericentre (optional). Used for break radius
- `quantiles::Vector{Real}`: The mass quantiles to calculate radii of
- `r_max::Real=1`: The maximum radius to calculate the velocity dispersion for
"""
function StellarScalars(snap::Snapshot;
        quantiles=[0.01, 0.03, 0.1, 0.5, 0.9, 0.99], 
        r_max=1,
        delta_t = NaN
    )

    # scalars
    bound_mass = sum(snap.weights[bound_particles(snap)])
    σ_v_1d = σv_1d(snap, r_max=r_max)
    r_break = break_radius(σ_v_1d, delta_t)
    r_quantile = find_r_quantile_star(snap, quantiles)

    return StellarScalars(
        delta_t = delta_t,
        r_break = r_break,
        sigma_v = σ_v_1d,
        r_max_sigma = r_max,
        bound_mass = bound_mass,
        r_quantiles = r_quantile,
        quantiles = quantiles,
        time = snap.time,
       )
end



"""
    scale(prof::StellarDensity3D, r_scale, m_scale)

Scale the profile by the given factors, returning a new
instance of StellarDensity3D. Note that the mass scale affects
only the stellar mass. 
"""
function scale(prof::StellarDensity3D, r_scale::Real, m_scale::Real=1)

    ρ_scale = m_scale / r_scale^3
    log_r_shift = log10(r_scale)

    return StellarDensity3D(
        log_r = prof.log_r .+ log_r_shift,
        log_r_bins = prof.log_r_bins .+ log_r_shift,
        rho = prof.rho * ρ_scale,
        counts = prof.counts,
        ess = prof.ess,
    )
end



"""
    scale(prof::StellarDensity3D, r_scale, m_scale, m_scale_pot)

Scale the profile by the given factors, returning a new
instance of StellarScalars. Note that the mass scale affects
only the stellar mass and the potential mass scale
affects the velocity dispersion.
"""
function scale(prof::StellarScalars, r_scale::Real,  m_scale::Real, m_scale_pot::Real)
    ρ_scale = m_scale / r_scale^3
    v_scale = sqrt(m_scale_pot / r_scale)

    return StellarScalars(
        delta_t = prof.delta_t,
        r_break = prof.r_break * r_scale, 
        sigma_v = prof.sigma_v * v_scale,
        r_max_sigma = prof.r_max_sigma * r_scale,
        bound_mass = prof.bound_mass * m_scale,
        r_quantiles = prof.r_quantiles .* r_scale,
        quantiles = prof.quantiles,
        time = prof.time,
       )

end



"""
    find_r_quantile_star(snap, p)

Find the radius(ii) containing a given fraction of the total stellar mass of the system
`p` may be a real or a vector of real representing the fractions.
"""
function find_r_quantile_star(snap::Snapshot, p)
    r = radii(snap)
    m = snap.weights ./ sum(snap.weights)

    return quantile(r, m, p)
end



"""
    σv_x(snap; r_max=1)

Calculate the velocity dispersion in the x-direction
for particles within r_max of the centre.
"""
function σv_x(snap; r_max=1)
    vs = x_velocity(snap)
    filt = radii(snap) .< r_max
    w = snap.weights[filt]
    vs = vs[filt]
    return std(vs, w)
end


"""
    σv_1d(snap; r_max=1)

Calculate the average 1D orthoganal velocity dispersion
for particles within r_max of the centre. 
Assums relative to snapshot center
"""
function σv_1d(snap; r_max=1)
    filt = radii(snap) .< r_max

    w = snap.weights[filt]
    v = speeds(snap)[filt]

    σ = sqrt( sum(v .^ 2 .* w) / sum(w) )

    return σ / √3
end



@doc raw"""
	break_radius(σ, delta_t)

Given a radial velocity dispersion σ and the time since last pericentre (all
code units) calculate the break radius in kpc (code units).

See peñarrubia et al. (200X) for the original equation.
```math
r_{\rm break} = C\,\sigma_{v}\,\Delta t
```
where $C=0.55$ is an empirical fit
"""
function break_radius(σ::Real, delta_t::Real; C::Real=0.55)
    return C * σ  * delta_t
end

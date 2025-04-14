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
    rho::Vector{F}

    "error in the density"
    rho_err::Vector{F}

    "number of particles in bin"
    counts::Vector{F}

    "effective sample size of bin (accounting for weights)"
    ess::Vector{F}

    annotations::Dict{String, Any} = Dict{String, Any}()
    "1D averaged velocity dispersion"
end


"""
    StellarDensity3D(snap; args...)

Calculate the 3D stellar density profile of the snapshot.

# Arguments
- `snap::Snapshot`: The snapshot to calculate the profile of
- `bins::Vector{Real}`: The bin edges (in log r) to use for the profile
- `quantiles::Vector{Real}`: The mass quantiles to calculate radii of
- `r_max::Real=1`: The maximum radius to calculate the velocity dispersion for
- `delta_t::Real=NaN`: The time since the last pericentre (optional). Used for break radius
"""
function StellarDensity3D(snap; delta_t=NaN, bins=nothing, 
        quantiles=[0.01, 0.03, 0.1, 0.5, 0.9, 0.99], 
        r_max=1,
    )

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

    # scalars
    σ_v_1d = σv_1d(snap, r_max=r_max)
    r_break = break_radius(σ_v_1d, delta_t)
    r_quantile = find_r_quantile_star(snap, quantiles)

    bound_mass = sum(snap.weights[bound_particles(snap)])

    return StellarDensity3D(
        log_r=log_r,
        log_r_bins=log_r_bins,
        rho=rho,
        rho_err=rho_err,
        counts = counts,
        ess = ess,
        annotations = Dict(
            "sigma_vx"=> σ_v_1d, 
            "quantiles"=> quantiles,
            "r_quantile"=> r_quantile,
            "time" =>  snap.time,
            "delta_t" => delta_t, 
            "r_break" => r_break,
            "bound_mass" => bound_mass,
           )
    )
end



"""
    scale(prof::StellarDensity3D, r_scale, v_scale, m_scale)

Scale the profile by the given factors, returning a new
instance of StellarDensity3D. Note that the mass scale affects
only the stellar mass. (Changes to the potential would change 
the v_scale but this only affects the velocity dispersion)
"""
function scale(prof::StellarDensity3D, r_scale::Real, v_scale::Real, m_scale::Real=1)

    ρ_scale = m_scale / r_scale^3
    log_r_shift = log10(r_scale)

    return StellarDensity3D(
        log_r = prof.log_r .+ log_r_shift,
        log_r_bins = prof.log_r_bins .+ log_r_shift,
        rho = prof.rho * ρ_scale,
        rho_err = prof.rho_err * ρ_scale,
        counts = prof.counts,
        ess = prof.ess,
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

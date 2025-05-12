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

    "snapshot time"
    time::F
end


break_radius(s::StellarScalars) = s.r_break
time_since_peri(s::StellarScalars) = s.delta_t
velocity_dispersion(s::StellarScalars) = s.sigma_v
mass(s::StellarScalars) = s.bound_mass


"""
    StellarScalars(snap::Snapshot)

Compute essential scalar properties of the snapshot.

- `delta_t::Real=NaN`: The time since the last pericentre (optional). Used for break radius
- `r_max::Real=1`: The maximum radius to calculate the velocity dispersion for
"""
function StellarScalars(snap::Snapshot;
        r_max=1,
        delta_t = NaN
    )

    # scalars
    bound_mass = sum(snap.weights[bound_particles(snap)])
    σ_v_1d = σv_1d(snap, r_max=r_max)
    r_break = break_radius(σ_v_1d, delta_t)

    return StellarScalars(
        delta_t = delta_t,
        r_break = r_break,
        sigma_v = σ_v_1d,
        r_max_sigma = r_max,
        bound_mass = bound_mass,
        time = snap.time,
       )
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
        time = prof.time,
       )

end

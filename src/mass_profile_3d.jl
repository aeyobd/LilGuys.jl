import TOML
import LinearAlgebra: ×, ⋅


# TODO: split MassProfile3D in half


"""
    MassProfile3D(<keyword arguments>)

A struct representing a 3-dimensional 
mass profile and cummulative / summative properties
(of a snapshot).

All properties are in code units.
"""
@kwdef struct MassProfile3D
    "Radius for v_circ"
    radii::Vector{F}

    "The number of particles within each radius"
    counts::Vector{F}

    "Circular velocity"
    v_circ::Vector{F}

    "Circular velocity error"
    v_circ_err::Vector{F}

    # quantiles (for convenience)
    "the radius at which q of the total mass is enclosed foreach quantile"
    r_quantile::Vector{F}

    "quantiles used for `r_quantile`"
    quantiles::Vector{F}

    # scalars
    """ Total Energy """
    E::F

    "Total Kinetic energy"
    K::F

    "Total Potential energy"
    W::F

    "Total Angular momentum vector wrt centre"
    L::Vector{F}

    "Maximum circular velocity"
    v_circ_max::F

    "Radius of maximum circular velocity"
    r_circ_max::F

    "Number of bound particles"
    N_bound::Int

    "Snapshot time"
    time::F = NaN
end



function Base.print(io::IO, prof::MassProfile3D)
    TOML.print(io, struct_to_dict(prof))
end


function MassProfile3D(filename::String)
    t = dict_to_tuple(TOML.parsefile(filename))
    return MassProfile3D(;t...)
end



"""
    MassProfile3D(snap; bins=100, filt_bound=true)

Given a snapshot, computes the density, circular velocity, and energy; 
returning an MassProfile3D object.

`bins` may be an integer, array, or function and is passed to `histogram`
"""
function MassProfile3D(snap::Snapshot;
        bins=nothing,
        filt_bound=true,
        quantiles=[0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999],
        skip = nothing,
    )

    if filt_bound
        filt = bound_particles(snap)
        snap = snap[filt]
    end

    N_bound = length(snap)
    if N_bound == 0
        throw(ArgumentError("No bound particles in snapshot"))
    end

    log_r_snap = log10.(radii(snap))
    if isnothing(skip)
        skip = min(round(Int, default_n_per_bin(log_r_snap) / 1.5), 200)
        @info "MassProfile3D binsize: $skip"
    end

    K = kinetic_energy(snap)

    if !isnothing(snap.potential)
        W = potential_energy(snap)
        E = W + K
    else
        E = NaN
        W = NaN
    end

    L = angular_momentum(snap)


    r_circ, v_c, n_circ = v_circ(snap, skip=skip)
    v_circ_err = v_c ./ sqrt.(n_circ)

    fit = fit_v_r_circ_max(r_circ, v_c)
    v_circ_max = fit.v_circ_max
    r_circ_max = fit.r_circ_max

    r_quantile = 10 .^ quantile(log_r_snap, quantiles)

    return MassProfile3D(
        radii=r_circ,
        counts=n_circ,
        v_circ=v_c,
        v_circ_err=v_circ_err,
        # quantiles
        r_quantile=r_quantile,
        quantiles=quantiles,
        # scalars
        v_circ_max=v_circ_max,
        r_circ_max=r_circ_max,
        E=E,
        K=K,
        W=W,
        L=L,
        N_bound=N_bound,
        time=snap.time,
    )
end


"""
    v_circ(snap; filter_bound=:recursive_1d, skip=10)

Returns a list of the sorted radii and circular velocity from a snapshot for the given centre.
Skips every `skip` particles (i.e. each radius contains n times skip particles inclusive)
filter_bound may be
- :recursive_1d: recursively removes unbound particles, recomputing the potential each time assuming spherical symmetry
- :recursive: like recursive_1d but computes the full potential
- :simple: only removes unbound particles once
- :false: does not filter particles
"""
function v_circ(snap::Snapshot; filter_bound=:recursive_1d, skip::Integer=10)
    r = radii(snap)
    m = snap.masses

    if filter_bound != :false
        if filter_bound == :simple
            filt = bound_particles(snap)
        elseif filter_bound == :recursive_1d
            filt = bound_particles_recursive_1D(snap)
        else
            throw(ArgumentError("Unknown filter_bound: $filter_bound"))
        end

        r = r[filt]
        m = m[filt]
    end

    if sum(r) == 0
        @warn "No bound particles in snapshot"
        return zeros(0), zeros(0), zeros(0)
    end

    idx = sortperm(r)
    m = snap.masses[idx]
    r = r[idx]
    M = cumsum(m)

    idx_skip = skip:skip:length(r)
    r = r[idx_skip]
    M = M[idx_skip]
    v_c = v_circ.(r, M)

    return r, v_c, idx_skip
end



"""
    fit_v_r_circ_max(r, v_circ; q=80, p0=[6., 30.])

Fits the maximum circular velocity of a rotation curve assuming a NFW
profile. Returns the parameters of the fit and the range of radii used.
The fit is done on the top `q` quantile of the circular velocity, as th
tails often extend strangely
"""
function fit_v_r_circ_max(r::AbstractArray{<:Real}, v_circ::AbstractArray{<:Real}; q=0.9, p0=nothing)
    if length(r) != length(v_circ)
        throw(DimensionMismatch("r and v_circ must have the same length. Got $(length(r)) and $(length(v_circ))"))
    elseif length(r) < 2
        @warn "r and v_circ must have at least 2 elements"
        return (; r_circ_max=NaN, v_circ_max=NaN, fit=nothing, converged=false, r_min=NaN, r_max=NaN)
    end

    if p0 == nothing
        idx = argmax(v_circ)
        p0 = [r[idx], v_circ[idx]]
    end


    filt = v_circ .> quantile(v_circ, q)

    local fit, converged

    if sum(filt) < 2
        converged = false
        fit = nothing
    else
        try
            fit = Interface.LsqFit.curve_fit(_v_circ_max_model, r[filt], v_circ[filt], p0)
            converged = fit.converged
        catch ArgumentError
            converged = false
            fit = nothing
        end
    end


    if !converged
        @warn "Fit did not converge, using simple maximum."
        idx = argmax(v_circ)
        v_circ_max = v_circ[idx]
        r_circ_max = r[idx]
        r_min = NaN
        r_max = NaN
    else
        r_circ_max = fit.param[1]
        v_circ_max = fit.param[2]
        r_min = minimum(r[filt]) 
        r_max = maximum(r[filt])
    end

    return (; 
        r_circ_max=r_circ_max,
        v_circ_max=v_circ_max,
        fit=fit,
        converged=converged,
        r_min=r_min,
        r_max=r_max
       )
end


"""
    fit_v_r_max(snap; kwargs...)

Fits circular velocity of snapshot
"""
function fit_v_r_circ_max(snap; kwargs...)
    r, v_c = v_circ(snap)
    return fit_v_r_circ_max(r, v_c; kwargs...)
end



@doc raw"""
NFW circular velocity model.

```math
v_circ(r) = v_circ_max/β * sqrt( (log(1 + x) - x/(1+x)) / x )
```

where x = r / r_circ_max * α_nfw
and α_nfw, β are the solution so v_circ_max_model(1) = 1.
"""
function _v_circ_max_model(r, param)
	r_circ_max, v_circ_max = param

	x = r ./ r_circ_max .* α_nfw
    β = 0.46499096281742197
    return @. v_circ_max / β  * sqrt(A_NFW(x) / x)
end


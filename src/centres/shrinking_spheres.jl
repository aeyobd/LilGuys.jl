import Base: @kwdef
import LinearAlgebra: norm


@kwdef struct _ShrinkingSpheresParams
    itermax::Int = 100
    r_cut_0::Float64 # see constructor
    r_factor::Float64 = 0.975
    N_min::Int # see constructor
    dx_atol::Float64 = 1e-2
    dx_rtol::Float64 = NaN
    dN_min::Int = 0
    x0::Vector{Float64} # see constructor
    mode::Symbol = :quantile
    r_min::Float64 = 0
    r_max::Float64 = Inf
    filter_unbound::Bool = true

    N::Int
end



"""
defines and initializes default parameters for the shrinking spheres algorithm
"""
function _ShrinkingSpheresParams(positions; f_min::Real=0.001, N_min::Real=100, 
        r_cut_0=nothing, x0=nothing, kwargs...)

    kwargs = Dict{Symbol,Any}(kwargs)

    N = size(positions, 2)
    kwargs[:N] = N


    if !isnothing(f_min)
        N_f_min = round(Int, f_min * N)
        if !isnothing(N_min)
            N_min = max(N_f_min, N_min)
        else
            N_min = N_f_min
        end
    end
    kwargs[:N_min] = N_min


    if isnothing(x0)
        x0 = centroid(positions)
    end
    kwargs[:x0] = x0

    if isnothing(r_cut_0)
        r_cut_0 = maximum(calc_r(positions, x0))
    end
    kwargs[:r_cut_0] = r_cut_0

    return _ShrinkingSpheresParams(;kwargs...)
end



@kwdef mutable struct SS_State <: AbstractState
    centre::Centre
    filt::BitVector
    kwargs::Dict{Symbol,Any}
end



"""
    SS_State(snap, kwargs...)

Given a snapshot, initialize a shrinking spheres state. 
Sets the initial centre to be from `centre_potential_percen` (i.e. using the 
centroid of the 5% most bound particles). 

See [`shrinking_spheres`](@ref) for arguments.

"""
function SS_State(snap::Snapshot; kwargs...)
    q0 = 0.05
    cen = centre_potential_percen(snap, q0)

    return SS_State(; centre=cen, filt=trues(size(snap.positions, 2)), kwargs=Dict{Symbol,Any}(kwargs))
end




"""
    calc_centre!(state::SS_State, snap::Snapshot)

Finds the centre using the shrinking sphere method
"""
function calc_centre!(state::SS_State, snap)
    # rearrange filter to match with snapshot's index (used for multiple iterations)
    idx = sortperm(snap.index)
    idx_r = invperm(idx)
    state.filt = state.filt[idx_r]

    params = _ShrinkingSpheresParams(snap.positions[:, state.filt]; state.kwargs...)


    _, filt = _shrinking_spheres(snap.positions[:, state.filt],
        params)

    state.filt[state.filt] .= filt

    state.centre = mean_centre(snap, state.filt)

    _filt_unbound!(state, snap, params)


    state.filt = state.filt[idx]
    return state
end



"""
    calc_next_centre!(state::SS_State, snap::Snapshot)

Finds the centre using the shrinking sphere method
"""
function calc_next_centre!(state::SS_State, snap::Snapshot)
    return calc_centre!(state, snap)
end



"""
    _bound_particles(state::SS_State, snap::Snapshot)

Retrieves a bitvector of particles that are bound with respect to the 
centre of the shrinking spheres state. 
Assumes that the potential of the snapshot is the real potential.
"""
function _bound_particles(state::SS_State, snap::Snapshot)
    v = calc_r(snap.velocities[:, state.filt] .- state.centre.velocity)

    Φs = snap.Φs[state.filt]
    ϵ = calc_E_spec.(Φs, v)
    filt_bound = ϵ .< 0
    
    return filt_bound
end


"""
    _filt_unbound!(state::SS_State, snap::Snapshot, params::_ShrinkingSpheresParams)

Filters out unbound particles from the state.
"""
function _filt_unbound!(state::SS_State, snap::Snapshot, params::_ShrinkingSpheresParams)
    if !params.filter_unbound
        return
    end
    filt_bound = _bound_particles(state, snap)
    if sum(filt_bound) < params.N_min
        @info "Too few particles are bound. "
        return
    end

    state.filt[state.filt] .= filt_bound

    Ncut = sum( @. !filt_bound)
    @info "Cut unbound particles: $Ncut" 

    state.centre = mean_centre(snap, state.filt)
end



"""
    shrinking_spheres(positions; <keyword arguments>)

Compute the centre using the iterative shrinking spheres method. Return the centre and a filter of the positions used to calculate the centre.

Implementation of the classic shrinking spheres method, described in Powers et
al. (2003, §2.5). The algorithm calculates the centroid and then removes
particles that are beyond a cutoff radius. The cutoff radius is decreased by a
fixed factor each step until the stopping criterion is met (reaching below a 
minimum number of particles or approximate convergence).
The final step removes unbound particles (if set) and recalculates the centre.

# Arguments
- `positions::AbstractMatrix{<:Real}`: The positions of the particles, with shape (3, N).
- `x0::AbstractVector{<:Real}`: The initial guess for the centroid. Defaults to the centroid of the positions using the state's filter
- `r_cut_0::Real`: The initial cutoff radius. Defaults to maximum radius of particles
- `r_factor::Real=0.975`: The factor by which the cutoff radius is decreased each iteration. `r_factor_` is interpreted as the quantile of particles to keep (if mode is quantile) or the ratio by which`r_max` decreases (if mode is ratio).
- `r_max::Real=Inf` Continue until the cutoff radius is below this value, ignoring other stopping criteria. Ignored if set to infinity.
- `r_min::Real=0`: Stop if the cutoff radius is below this value.
- `f_min::Real=0.001`: Stop if the number of particles is below f_min * N.
- `N_min::Int=100`: Stop if the number of particles is below N_min.
- `itermax::Int=100`: Maximum number of iterations.
- `dx_atol::Real=1e-6`: Stop if the absolute change in the centroid is below this value.
- `dx_rtol::Real=NaN`: Stop if the relative change in the centroid is below this value.
- `dN_min::Int=0`: Stop if the number of particles removed in a step is below this value.
- `mode::Symbol=:quantile`: The mode for determining the cutoff radius. Either :quantile or :ratio.
- `filter_unbound::Bool=true`: If true, remove unbound particles after each iteration. Uses precomputed potentials.
"""
function shrinking_spheres(positions::AbstractMatrix{T}; kwargs...) where T<:Real
    params = _ShrinkingSpheresParams(positions; kwargs...)
    _check_params(positions, params)
    return _shrinking_spheres(positions, params)
end




function _shrinking_spheres(positions, params::_ShrinkingSpheresParams)
    r_cut = params.r_cut_0
    x0 = copy(params.x0)
    filt = trues(params.N)

    status = nothing

    for i in 1:params.itermax
        radii = calc_r(positions[:, filt], x0)

        if params.mode == :quantile
            r_cut = quantile(radii, params.r_factor)
        else
            r_cut *= params.r_factor
        end

        filt_r = radii .< r_cut

        xnew = centroid(positions[:, filt][:, filt_r])

        dx = calc_r(xnew, x0)
        N = sum(filt)
        dN = sum( @. !filt_r)
        

        _print_status(i, x0, N, dx, dN, r_cut)

        status = _is_complete(params, x0, N, dx, dN, r_cut)

        if !isnothing(status)
            break
        end

        # only update if we continue
        x0 = xnew
        filt[filt] .= filt_r
    end

    if status == nothing
        status = "itermax"
    end

    @info "Shrinking spheres completed with status $status"

    return x0, filt
end


function _check_params(positions, params::_ShrinkingSpheresParams)
    if size(positions, 2) != params.N
        throw(ArgumentError("Number of particles in positions does not match N"))
    elseif params.r_cut_0 < 0
        throw(ArgumentError("r_cut_0 must be positive"))
    elseif params.r_factor <= 0 || params.r_factor >= 1
        throw(ArgumentError("r_factor must be in (0, 1)"))
    elseif params.itermax <= 0
        throw(ArgumentError("itermax must be positive"))
    elseif params.mode ∉ [:quantile, :ratio]
        throw(ArgumentError("mode must be :quantile or :ratio"))
    end
end


function _is_complete(params::_ShrinkingSpheresParams, x0, N, dx, dN, r_cut)
    if N < params.N_min
        return "N"
    end

    if params.r_max < Inf 
        if r_cut < params.r_max
            return "r_max"
        else
            return nothing
        end
    end

    if dN < params.dN_min
        return "dN"
    elseif dx < params.dx_atol && (dN > 0)
        return "dx"
    elseif (!isnan(params.dx_rtol)
            && dx / norm(x0) < params.dx_rtol)
        return "dx_rel"
    elseif r_cut < params.r_min
        return "r_min"
    else
        return nothing
    end
end


function _print_status(i, x0, N, dx, dN, r_cut)
    @info "i=$i, $x0, N=$N, dx=$dx, dN=$dN, r_cut=$r_cut"
end


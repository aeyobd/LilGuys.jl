import Base: @kwdef
import StatsBase: percentile


@kwdef mutable struct SS_State <: AbstractState
    centre::Centre
    f_min::Float64 = 0.1
    percen_i::Float64 = 95
    percen::Float64 = 99
    itermax::Int = 100

    filt::BitVector

    iter::Int = 0
    verbose::Bool = false
    Np::Int = 0
end



"""
    SS_State(snap, kwargs...)

Given a snapshot, initialize a shrinking spheres state. 

Parameres
---------

"""
function SS_State(snap::Snapshot; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    kwargs[:filt] = trues(length(snap))

    cen = centre_potential_percen(snap)
    kwargs[:centre] = cen
    kwargs[:Np] = length(snap)

    return SS_State(;kwargs...)
end



function Base.copy(state::SS_State)
    kwargs = Dict{Symbol,Any}()
    for k in fieldnames(SS_State)
        v = getfield(state, k)
        if v !== nothing
            kwargs[k] = copy(v)
        end
    end
    return SS_State(;kwargs...)
end



function calc_next_centre!(state::SS_State, snap::Snapshot)
    filt = sortperm(snap.index)
    snap = snap[filt]
    _shrink_sphere!(state, snap)
    if state.verbose
        println("done")
    end
    return state
end



"""
Finds the centre using the shrinking sphere method
"""
function calc_centre!(state::SS_State, snap)
    _print_header(state)
    filt = sortperm(snap.index)
    snap = snap[filt]
    state.iter = 0

    for i in 1:state.itermax
        dN = _shrink_sphere!(state, snap)
        state.iter += 1
        if dN == 0
            break
        end
    end

    if state.verbose
        println("completed")
    end
    return state
end




function is_done(state::SS_State)
    if state.iter > state.itermax
        println("Warning: max iterations reached")
        return true
    end
    N_min = length(state.filt) * state.f_min
    return sum(state.filt) <= N_min
end




function _shrink_sphere!(state::SS_State, snap::Snapshot)
    filt_i = state.filt
    rs = calc_r(snap.positions[:, filt_i] .- state.centre.position)
    v = calc_r(snap.velocities[:, filt_i] .- state.centre.velocity)

    Φs = calc_radial_discrete_Φ(snap.masses[filt_i], rs)
    ϵ = calc_E_spec.(Φs, v)
    filt_bound = ϵ .< 0

    r_cut = percentile(rs, state.percen)
    filt_r = rs .<= r_cut

    dN_bound = _add_filter!(state, filt_bound)
    if dN_bound == 0
        filt_bound = trues(length(filt_bound))
    end
    filt_r = filt_r[filt_bound]

    dN_r = _add_filter!(state, filt_r)

    if state.verbose
        println("$dN_bound,$dN_r,$r_cut")
    end
    masses = snap.masses[state.filt]
    state.centre = weighted_centre(snap[state.filt], masses)

    return dN_r + dN_bound
end

function _print_header(state)
    if state.verbose
        println("dN_bound,dN_r,r_cut")
    end
end



function _add_filter!(state, filt)
    N  = sum(state.filt)
    N_after = sum(filt)
    if N_after <= state.f_min * state.Np
        return 0
    end
    state.filt[state.filt] .= filt

    return N - N_after
end


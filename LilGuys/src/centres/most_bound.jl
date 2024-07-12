import Base: @kwdef
import StatsBase: percentile

not = !


@kwdef mutable struct MostBoundState <: AbstractState
    centre::Centre
    f_initial::Float64 = 0.01
    f_min::Float64 = 0.001

    itermax_initial::Int = 20
    itermax::Int = 5
    percen::Float64 = 90

    r_max::Float64 = 3.0

    filt::BitVector

    verbose::Bool = false

    Np::Int
    i::Int = 1
end




"""
    MostBoundState(snap, kwargs...)

Given a snapshot, initialize a state

"""
function MostBoundState(snap::Snapshot; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    kwargs[:filt] = trues(length(snap))
    kwargs[:centre] = centre_potential_percen(snap, 5)
    kwargs[:Np] = length(snap)
    println("Np = ", kwargs[:Np])

    return MostBoundState(;kwargs...)
end




function calc_next_centre!(state::MostBoundState, snap::Snapshot)
    idx = sortperm(snap.index)
    idx_r = invperm(idx)
    state.i += 1

    # filter in snapshot ordering
    state.filt = state.filt[idx_r]
    state.centre = weighted_centre(snap[state.filt], snap.masses[state.filt])

    for i in 1:state.itermax
        cont =  _cut_unbound!(state, snap)
        if !cont
            break
        end
    end

    if state.verbose
        println("completed $(state.i)th centre")
        println("number of particles: ", sum(state.filt))
    end

    # return to original ordering
    state.filt = state.filt[idx]
    return state
end



"""
Finds the centre using the shrinking sphere method
"""
function calc_centre!(state::MostBoundState, snap)
    _print_header(state)

    if state.verbose
        println("using initial centre ", state.centre)
    end

    for i in 1:state.itermax_initial
        cont = _shrink_sphere!(state, snap)
        if !cont
            break
        end
    end

    # return to original ordering
    state.filt = state.filt[sortperm(snap.index)]

    if state.verbose
        println("completed")
        println("number of particles: ", sum(state.filt))
    end

    return state
end


function _cut_unbound!(state::MostBoundState, snap::Snapshot)
    v = calc_r(snap.velocities[:, state.filt] .- state.centre.velocity)
    ϵ = calc_E_spec.(snap.Φs[state.filt], v)
    ϵ_cut = 0

    filt = ϵ .< ϵ_cut
    r = calc_r(snap.positions[:, state.filt] .- state.centre.position)
    filt .&= r .< state.r_max

    state.filt[state.filt] .= filt

    println("dN_bound = ", sum(not.(filt)))

    state.centre = weighted_centre(snap[state.filt], snap.masses[state.filt])
    return sum(not.(filt)) > 0
end



function _shrink_sphere!(state::MostBoundState, snap::Snapshot)
    filt_i = state.filt

    v = calc_r(snap.velocities[:, filt_i] .- state.centre.velocity)
    Φs = snap.Φs[filt_i]
    ϵ = calc_E_spec.(Φs, v)

    ϵ_cut = percentile(ϵ, state.percen)

    filt = ϵ .< ϵ_cut

    dN = sum(not.(filt))

    dN = sum(not.(filt))

    if state.verbose
        println("dN = $dN")
    end

    state.filt[filt_i] .= filt

    state.centre = weighted_centre(snap[state.filt], snap.masses[state.filt])

    println("Nmax ", state.f_initial * state.Np)
    return sum(state.filt) > state.f_initial * state.Np
end



function _print_header(state::MostBoundState)
    if state.verbose
        println("dN_bound,dN_r,dN")
    end
end

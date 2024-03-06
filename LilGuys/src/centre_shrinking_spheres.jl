import Base: @kwdef
import StatsBase: percentile


@kwdef mutable struct SS_State
    f_min::Float64 = 0.1
    const percen::Float64 = 95

    filt::BitVector
    Φs::Vector{Float64}
    rs::Vector{Float64}

    x_c::Vector{Float64}
    v_c::Vector{Float64}
    iter::Int = 0
    itermax::Int = 100

    history = nothing
end


@kwdef mutable struct SS_History
    xs::Vector{Vector{Float64}} = []
    vs::Vector{Vector{Float64}} = []
    Rmax::Vector{Float64} = []
    N::Vector{Int} = []
    dN_R::Vector{Int} = []
    dN_bound::Vector{Int} = []
end


function SS_State(snap::Snapshot; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)

    kwargs[:x_c] = zeros(3)
    kwargs[:v_c] = zeros(3)
    kwargs[:Φs] = zeros(length(snap))
    kwargs[:rs] = calc_r(snap.positions)
    kwargs[:filt] = trues(length(snap))

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




function ss_centre(out::Output; f_step_factor=1, kwargs...)
    snap = out[1]
    filt = sortperm(snap.index)
    state = ss_centre(snap[filt]; kwargs...)

    states = [state]

    xs = Matrix{Float64}(undef, length(out), 3)
    vs = Matrix{Float64}(undef, length(out), 3)
    r_cut = 10

    for snap in out[2:end]
        filt = sortperm(snap.index)
        snap = snap[filt]
        state = copy(state)

        filter_unbound!(state, snap)
        filter_radius!(state, snap, r_cut)
        filter_unbound!(state, snap)

        push!(states, state)
    end

    return states
end




"""
Finds the centre using the shrinking sphere method
"""
function ss_centre(snap::Snapshot; verbose=false, history=false, kwargs...)
    state = SS_State(snap; kwargs...)
    update_centre!(state, snap)

    while !is_done(state)
        dN_r, r_cut = filter_edges!(state, snap)
        update_centre!(state, snap)
        dN_bound = filter_unbound!(state, snap)

        update_centre!(state, snap)
        if verbose
            println("iteration ", state.iter)
            println("r_cut = ", r_cut)
            println("dN_r = ", dN_r)
            println("dN_bound = ", dN_bound)
            println("N = ", sum(state.filt))
        end
        if history
            update_history!(state, r_cut, dN_r, dN_bound)
        end

        if dN_r == 0 && dN_bound == 0
            break
        end
        state.iter += 1
    end

    return state
end


function update_history!(state::SS_State, r_cut, dN_r, dN_bound)
    if state.history == nothing
        state.history = SS_History()
    end
    push!(state.history.N, sum(state.filt))
    push!(state.history.xs, state.x_c)
    push!(state.history.vs, state.v_c)
end



function is_done(state::SS_State)
    if state.iter > state.itermax
        println("Warning: max iterations reached")
        return true
    end
    N_min = length(state.filt) * state.f_min
    return sum(state.filt) <= N_min
end


function filter_edges!(state::SS_State, snap)
    r_cut = percentile(state.rs, state.percen)

    dN = filter_radius!(state, snap, r_cut)

    if state.history !== nothing
        push!(state.history.Rmax, r_cut)
        push!(state.history.dN_R, dN)
    end
    return dN, r_cut
end


"""
Filters particles beyond the specified radius
"""
function filter_radius!(state::SS_State, snap::Snapshot, r_cut)
    filt_i = state.filt
    filt = state.rs .<= r_cut
    N = sum(filt_i)
    N_after = sum(filt)

    if N_after < N * state.f_min
        return 
    end
    state.filt[filt_i] .= filt
    update_centre!(state, snap)

    return N - N_after
end


function filter_unbound!(state::SS_State, snap::Snapshot)
    filt_i = state.filt
    N  = sum(filt_i)

    r = state.rs
    v = calc_r(snap.velocities[:, filt_i] .- state.v_c)

    recalc_Φ!(state, snap)
    ϵ = calc_E_spec.(state.Φs, v)

    filt = ϵ .< 0

    N_after = sum(filt)

    if N_after < N * state.f_min
        update_centre!(state, snap)
        return 
    end

    state.filt[filt_i] .= filt

    if state.history !== nothing
        dN_bound = N - N_after
        push!(state.history.dN_bound, dN_bound)
    end

    update_centre!(state, snap)

    return N - N_after
end


function recalc_Φ!(state::SS_State, snap)
    ms = snap.masses[state.filt]
    rs = state.rs
    state.Φs = calc_radial_discrete_Φ(ms, rs)
end


function update_centre!(state::SS_State, snap::Snapshot)
    state.x_c = centroid(snap.positions[:, state.filt])
    state.v_c = centroid(snap.velocities[:, state.filt])
    state.rs = calc_r(snap.positions[:, state.filt] .- state.x_c)
end



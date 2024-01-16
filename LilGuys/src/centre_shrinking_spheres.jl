import Base: @kwdef
import StatsBase: percentile


@kwdef mutable struct SS_State
    const N_min::Int
    const percen::Float64 = 95
    const snap::Snapshot # doesn't enforce no changes....

    Φ_func = calc_radial_Φ(snap.positions, snap.masses)

    filt::BitVector = trues(length(snap))
    x_c::Vector{Float64}
    v_c::Vector{Float64}
    iter::Int = 0

    # running info, could probably move to a separate struct and be optional
    # (but im grad student loll)
    xs::Vector{Vector{Float64}} = []
    vs::Vector{Vector{Float64}} = []
    Rmax::Vector{Float64} = []
    N::Vector{Int} = []
    dN_R::Vector{Int} = []
    dN_bound::Vector{Int} = []
end


function SS_State(snap::Snapshot; kwargs...)
    kwargs = Dict{Symbol,Any}(kwargs)
    if !haskey(kwargs, :f_min)
        f_min = 0.01
    else
        fmap = kwargs[:f_min]
        delete!(kwargs, :f_min)
    end
    kwargs[:N_min] = ceil(f_min * length(snap))

    kwargs[:snap] = snap
    x_c = centroid(snap.positions)
    v_c = centroid(snap.velocities)
    kwargs[:x_c] = x_c
    kwargs[:v_c] = v_c

    return SS_State(;kwargs...)
end



"""
Finds the centre using the shrinking sphere method
"""
function ss_centre(snap::Snapshot; verbose=false, kwargs...)
    state = SS_State(snap; kwargs...)

    while !is_done(state)
        filter_edges!(state)
        recalc_Φ!(state)
        filter_unbound!(state)
        update_centre!(state)
        if verbose
            println("iteration ", state.iter)
        end
    end

    return state
end


function update_centre!(state::SS_State)
    state.iter += 1

    push!(state.xs, state.x_c)
    push!(state.vs, state.v_c)
    state.x_c = centroid(state.snap.positions[:, state.filt])
    state.v_c = centroid(state.snap.velocities[:, state.filt])

    push!(state.N, sum(state.filt))
end


function is_done(state::SS_State)
    return sum(state.filt) <= state.N_min
end


function filter_edges!(state::SS_State)
    filt_i = state.filt
    r = calc_r(state.snap.positions[:, filt_i] .- state.x_c)
    r_cut = percentile(r, state.percen)

    filt = r .<= r_cut
    N = sum(filt_i)
    N_after = sum(filt)
    push!(state.Rmax, r_cut)
    push!(state.dN_R, N - N_after)
    state.filt[filt_i] .= filt
end


function filter_unbound!(state::SS_State)
    filt_i = state.filt
    N  = sum(filt_i)

    r = calc_r(state.snap.positions[:, filt_i] .- state.x_c)
    v = calc_r(state.snap.velocities[:, filt_i] .- state.v_c)

    ϵ = calc_E_spec.(state.Φ_func.(r), v)

    filt = ϵ .< 0
    N_after = sum(filt)
    if N_after < state.N_min
        filt = trues(length(filt))
        N_after = N
    end

    state.filt[filt_i] .= filt

    dN = N - N_after
    push!(state.dN_bound, dN)
end


function recalc_Φ!(state::SS_State)
    r_vecs = state.snap.positions .- state.x_c
    ms = state.snap.masses[state.filt]
    r_vecs = r_vecs[:, state.filt]
    state.Φ_func = calc_radial_Φ(r_vecs, ms)
end

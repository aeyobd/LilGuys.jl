import Base: @kwdef
import StatsBase: percentile

@kwdef mutable struct SS_State
    const N_min::Int
    const percen::Float64 = 95

    snap::Snapshot
    x_c::Vector{Float64}
    v_c::Vector{Float64}
    iter::Int = 1

    # running info
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

    kwargs[:snap] = copy(snap)
    x_c = centroid(snap.pos)
    v_c = centroid(snap.vel)
    kwargs[:x_c] = x_c
    kwargs[:v_c] = v_c

    return SS_State(;kwargs...)
end



"""
Finds the centre using the shrinking sphere method
"""
function ss_centre(snap::Snapshot; kwargs...)
    state = SS_State(snap; kwargs...)

    while !is_done(state)
        cut_outside!(state)
        recalc_Φ!(state)
        cut_unbound!(state)
        update_centre!(state)
    end

    return state
end


function update_centre!(state::SS_State)
    state.iter += 1

    push!(state.xs, state.x_c)
    push!(state.vs, state.v_c)
    state.x_c = centroid(state.snap.pos)
    state.v_c = centroid(state.snap.vel)

    push!(state.N, length(state.snap))
end


function is_done(state::SS_State)
    return length(state.snap) <= state.N_min
end


function cut_outside!(state::SS_State)
    r = calc_r(state.snap.pos .- state.x_c)
    r_cut = percentile(r, state.percen)

    filt = r .<= r_cut
    state.snap = state.snap[filt]

    N_cut = sum(map(!, filt))
    push!(state.Rmax, r_cut)
    push!(state.dN_R, length(state.snap))
end


function cut_unbound!(state::SS_State)
    x1 = state.snap.pos .- state.x_c
    v1 = state.snap.vel .- state.v_c

    E_kin = 0.5 * calc_r(v1).^2
    E_spec = state.snap.Φ .+ E_kin

    filt = E_spec .< 0
    N_cut = sum(map(!, filt))

    state.snap = state.snap[filt]
    push!(state.dN_bound, N_cut)
end


function recalc_Φ!(state::SS_State)
    pos = state.snap.pos .- state.x_c
    state.snap.Φ .= calc_radial_Φ(pos, state.snap.m).(calc_r(pos))
end

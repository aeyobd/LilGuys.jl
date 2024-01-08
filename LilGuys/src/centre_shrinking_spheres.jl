import Base: @kwdef

@kwdef mutable struct SS_State
    snap::Snapshot
    x_c::Vector{Float64}
    v_c::Vector{Float64}
    const N_min::Int
    δxs::Vector{Vector{Float64}} = []
    δvs::Vector{Vector{Float64}} = []
    iter::Int = 1
end


function SS_State(snap::Snapshot, f_min=0.2)
    kwargs = Dict{Symbol, Any}()
    kwargs[:snap] = copy(snap)
    kwargs[:N_min] = ceil(f_min * length(snap))
    x_c, v_c = centroid(snap.pos)
    kwargs[:x_c] = x_c
    kwargs[:v_c] = v_c

    return SS_State(;kwargs...)
end


function update!(state::SS_State, δx, δv)
    state.iter += 1
    push!(state.δxs, δx)
    push!(state.δvs, δv)

    state.x_c += state.δxs[end]
    state.v_c += state.δvs[end]

    state.snap.pos .-= δx
    state.snap.vel .-= δv
end


function is_done(state::SS_State)
    return length(state.snap) <= state.N_min
end

"""
Finds the centre using the shrinking sphere method
"""
function ss_centre(snap::Snapshot, perc=95, f_min=0.2)
    state = SS_State(snap, f_min)

    while !is_done(state)
        cut_outside!(state.snap, state.x_c, perc)
        calc_Φ!(state.snap)
        cut_unbound!(state.snap, state.x_c, state.v_c)

        update!(state, perc)
    end

    return state.x_c, state.v_c, state
end



function cut_outside(snap, cerc=95)
    r = r(snap.pos)
    r_cut = percentile(r, perc)
    println(r_cut)
    filt = r .<= r_cut
    return snap[filt]
end


function cut_unbound!(snap::Snapshot)
    x1 = snap.pos .- x0
    v1 = snap.vel .- v0

    E_kin = 0.5 * r(v1).^2
    E_spec = snap.Φ .+ E_kin

    filt = E_spec .< 0
    return snap[filt]
end


function potential_centre(snap::AbstractSnapshot; percen=5)
    threshhold = percentile(snap.Φ, percen)
    filt = snap.Φ .< threshhold
    return centroid(snap[filt])
end


function centroid(snap::AbstractSnapshot)
    pos_c, δr = centroid(snap.pos)
    vel_c, δv = centroid(snap.vel)

    return FuzzyPhase(pos_c, vel_c, δr, δv)
end

function centroid(snap::AbstractSnapshot, weights::Vector{F})
    pos_c = centroid(snap.pos, weights)
    vel_c = centroid(snap.vel, weights)
    δr = centroid_err(snap.pos .- pos_c, weights)
    δv = centroid_err(snap.vel .- vel_c, weights)

    return FuzzyPhase(pos_c, vel_c, δr, δv)
end


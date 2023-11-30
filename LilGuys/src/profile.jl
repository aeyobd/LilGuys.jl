# represents the 1D profile of some star
struct Profile
    r
    m
    M
    ρ_s
    ρ_dm
    ψ
end

function Profile(snap::Snapshot, r_s, n)
    r = norm.(snap.pos, 2)
    m = snap.m
    M = cumsum(m)
    ρ_s = ρ_s(r, r_s, n)
    ρ_dm = ρ_dm(r, r_s, n)
    ψ = ψ(r, r_s, n)
    return Profile(r, m, M, ρ_s, ρ_dm, ψ)
end


function get_r(x::Matrix{F})
    return reshape(sqrt.(sum(x.^2, dims=1)), :)
end

struct NFWParams
    x0
    r_s
end

function p_nfw(x::Matrix{F}, params::NFWParams)
    r = get_r(x .- params.x0)
    return @. 1/params.r_s^3 / (r/params.r_s) / (1 + r/params.r_s)^2
end

function get_center(snap::Snapshot)
    idx = argmin(snap.Φ)
    return snap.pos[:, idx]
end

function center(snap::Snapshot)
    snap1 = copy(snap)
    snap1.pos .-= get_center(snap)
    return snap1
end

function log_likelyhood_nfw(pos, params::NFWParams)
    if params.r_s < 0
        return Inf
    end
    return -sum(log.(p_nfw(pos, params)))
end

get_r(snapshot::Snapshot) = get_r(snapshot.pos)

# check these!
ρ_s(r, r_s, n) = exp(-(r/r_s)^n)
ρ_s_int(r, r_s, n) = -r_s^n * exp(-(r/r_s)^n) * (r/r_s)^n * (n+1)

function ρ_E21(x, x0, r_s, r_t, c)
end

function get_bound(snap::Snapshot)
    return snap[E_local(snap) .< 0]
end

V_circ(M, r) = r > 0 ? sqrt(G*M/r) : 0


function stellar_probabilities(snap::Snapshot, r_s, n)
end

function V_circ(snap::Snapshot, r)
    rs = get_r(snap)
    M = sum(rs < r) * snap.m
    return V_circ(M, r)
end

function most_bound(snap::Snapshot, percentile=0.2)

end

function KE(snap::Snapshot, v0=zeros(3))
    v0 = reshape(v0, 3, 1)
    return 0.5 * snap.m .* get_r(snap.vel .- v0).^2
end

function E_local(snap::Snapshot, v0=zeros(3))
    return KE(snap, v0) .+ snap.m .* snap.Φ
end

function E_tot(snap::Snapshot, v0=zeros(3))
    return sum(KE(snap, v0) .+ 0.5*snap.m .* snap.Φ .+ snap.m .* snap.Φ_ext)
end


function V_circ_max(snap::Snapshot)
    r = get_r(snap.pos)
    M = arange(length(snapshot)) .* snapshot.m
    return max(V_circ(M, r))
end


function angular_momentum(snap::Snapshot)
    return cross(snap.pos, snap.vel) .* snap.m
end


function sort_by_r(snap::Snapshot)
    return snap[sortperm(get_r(snap.pos))]
end


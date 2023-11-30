import StatsBase: percentile, mean, std


function find_centre(snap::Snapshot, perc=95, min_frac=0.2)
    N = length(snap)
    N_min = ceil(min_frac * N)

    snap1 = copy(snap)
    idx0 = argmin(snap1.Φ)
    xs = []
    vs = []


    x0 = snap1.pos[:, idx0]
    v0 = snap1.vel[:, idx0]

    while length(snap1) > N_min
        push!(xs, x0)
        push!(vs, v0)

        snap1 = cut_outside(snap1, x0, perc)
        calc_Φ!(snap1)
        snap1 = cut_unbound(snap1, x0, v0)

        x0 = centre_of(snap1.pos)
        v0 = centre_of(snap1.vel)
    end

    print(std(xs))
    print(std(vs))

    return x0, v0
end

function centre_of(x::Matrix{F})
    return mean(x, dims=2)
end


function calc_Φ!(snap::Snapshot)
    for i in 1:length(snap)
        r = get_r(snap.pos[:, 1:end .!= i] .- snap.pos[i])
        snap.Φ[i] = -snap.m * sum(1 ./ r)
    end
    return snap
end


function cut_outside(snap, c, perc=95)
    x1 = snap.pos .- c
    r = get_r(x1)
    r_cut = percentile(r, perc)
    println(r_cut)
    filt = r .<= r_cut
    return snap[filt]
end


function cut_unbound(snap::Snapshot, x0, v0)
    x1 = snap.pos .- x0
    v1 = snap.vel .- v0

    E_kin = 0.5 * get_r(v1).^2
    E_spec = snap.Φ .+ E_kin

    filt = E_spec .< 0
    return snap[filt]
end

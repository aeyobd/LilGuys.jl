

"""
Finds the centre using the shrinking sphere method
"""
function ss_centre(snap::Snapshot, perc=95, min_frac=0.2)
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

        x0, v0 = centroid(snap1.pos)
    end

    return x0, v0
end



function cut_outside(snap, cerc=95)
    r = r(snap.pos)
    r_cut = percentile(r, perc)
    println(r_cut)
    filt = r .<= r_cut
    return snap[filt]
end


function cut_unbound(snap::Snapshot, x0, v0)
    x1 = snap.pos .- x0
    v1 = snap.vel .- v0

    E_kin = 0.5 * r(v1).^2
    E_spec = snap.Φ .+ E_kin

    filt = E_spec .< 0
    return snap[filt]
end



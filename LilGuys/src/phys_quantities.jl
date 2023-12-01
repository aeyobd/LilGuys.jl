
function get_r(x::Matrix{F})
    return reshape(sqrt.(sum(x.^2, dims=1)), :)
end


function get_bound(snap::Snapshot)
    return snap[E_local(snap) .< 0]
end

V_circ(M, r) = r > 0 ? sqrt(G*M/r) : 0


function E_spec_kin(snap::Snapshot, v0=zeros(3))
    v0 = reshape(v0, 3, 1)
    return 0.5 .* get_r(snap.vel .- v0).^2
end

function E_spec(snap::Snapshot)
    return E_spec_kin(snap) .+ snap.Φ
end

function E_tot(snap::Snapshot, v0=zeros(3))
    return sum(snap.m .* E_spec_kin(snap, v0) .+ 0.5*snap.m .* snap.Φ .+ snap.m .* snap.Φ_ext)
end


function angular_momentum(snap::Snapshot)
    return cross(snap.pos, snap.vel) .* snap.m
end


using LinearAlgebra: ×

"""
The magnitude of a matrix of 3-vectors of shape 3×N
"""
function r(x::Matrix{F})
    return reshape(sqrt.(sum(x.^2, dims=1)), :)
end


function r(p::Particle) 
    return r(p.pos)
end


function v(p::Particle) 
    return r(p.vel)
end


function r(x::Vector{F})
    if length(x) != 3
        error("Vector must have length 3")
    end
    return sqrt(sum(x.^2))
end


function E_spec_kin(snap::Snapshot, v0=zeros(3))
    v0 = reshape(v0, 3, 1)
    return 0.5 .* r(snap.vel .- v0).^2
end


"""
Specific energy of each particle
"""
function E_spec(snap::Snapshot)
    return E_spec_kin(snap) .+ snap.Φ
end


function E_spec(Φ::F, v::F)
    return 0.5v^2 .+ Φ
end


function E_tot(snap::Snapshot, v0=zeros(3))
    return sum(snap.m .* E_spec_kin(snap, v0) .+ 0.5*snap.m .* snap.Φ .+ snap.m .* snap.Φ_ext)
end


function angular_momentum(snap::Snapshot)
    L = Matrix{F}(undef, 3, length(snap))
    for i in 1:length(snap)
        p = snap[i]
        L[:, i] .= p.pos × p.vel .* p.m
    end
    return L
end


"""
The circular velocity at radius r from the center of a mass M.
"""
function V_circ(M::F, r::F)
    if r == 0
        return 0
    elseif r < 0
        error("r must be positive")
    end
    return  sqrt(G*M/r)
end


function get_bound(snap::Snapshot)
    return snap[E_spec(snap) .< 0]
end

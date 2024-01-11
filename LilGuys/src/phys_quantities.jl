using LinearAlgebra: ×

"""
The magnitude of a matrix of 3-vectors of shape 3×N
"""
function calc_r(x::Matrix{T}) where T<:Real
    return reshape(sqrt.(sum(x.^2, dims=1)), :)
end



function calc_r(x::Vector{T}) where T<:Real
    if length(x) != 3
        error("Vector must have length 3")
    end
    return sqrt(sum(x.^2))
end


function calc_E_spec_kin(snap::Snapshot, v0=zeros(3))
    v0 = reshape(v0, 3, 1)
    r = calc_r(snap.vel .- v0)
    return 0.5 .* r.^2
end


"""
Specific energy of each particle
"""
function calc_E_spec(snap::Snapshot)
    return calc_E_spec_kin(snap) .+ snap.Φ
end


function calc_E_spec(Φ::Real, v::Real)
    return 0.5v^2 .+ Φ
end


function calc_E_tot(snap::Snapshot, v_vec_0=zeros(3))
    ms = snap.masses
    return sum(ms .* E_spec_kin(snap, v_vec) 
               .+ 1/2*ms .* snap.Φs 
               .+ ms .* snap.Φs_ext)
end


function calc_L(snap::Snapshot)
    L = Matrix{F}(undef, 3, length(snap))

    for i in 1:length(snap)
        p = snap[i]
        x_vec = p.positions[:, i]
        v_vec = p.velocities[:, i]
        m = p.masses[i]
        L[:, i] .= m .* (x_vec × v_vec)
    end
    return L
end

function calc_L_tot(snap::Snapshot)
    return sum(calc_L(snap), dims=2)
end


"""
The circular velocity at radius r from the center of a mass M.
"""
function calc_V_circ(M::Real, r::Real)
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

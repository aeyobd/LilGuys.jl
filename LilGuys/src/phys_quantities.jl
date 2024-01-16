using LinearAlgebra: ×
using EllipsisNotation

"""
The magnitude of a matrix of 3-vectors of shape 3×N
"""
function calc_r(x::Matrix{T}) where T<:Real
    if size(x, 1) != 3
        error("Matrix must have shape 3×N")
    end
    return _calc_r(x)
end


function calc_r(a::Array{T}, b::Array{T}) where T<:Real
    return calc_r(a .- b)
end

function calc_r(x::Vector{T}) where T<:Real
    if length(x) != 3
        error("Vector must have length 3")
    end
    #return sqrt(sum(x.^2))
    return _calc_r(x)
end


function _calc_r(x::Array{T}) where T<:Real
    r = sqrt.(sum(x.^2, dims=1))
    return r[1, ..] #sum doesn't actually reduce the dimension
end


function calc_r(snap::Snapshot)
    return calc_r(snap.positions)
end


function calc_v(snap::Snapshot)
    return calc_r(snap.velocities)
end

function calc_v(snap::Snapshot, v0::Vector{T}) where T<:Real
    return calc_r(snap.velocities .- v0)
end

function calc_E_spec_kin(snap::Snapshot, v0=zeros(3))
    v = calc_v(snap, v0)
    return 0.5 .* v.^2
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


function calc_L_spec(snap::Snapshot)
    L = Matrix{F}(undef, 3, length(snap))

    for i in 1:length(snap)
        x_vec = snap.positions[:, i]
        v_vec = snap.velocities[:, i]
        L[:, i] .= x_vec × v_vec
    end
    return L
end

function calc_L(snap::Snapshot)
    return calc_L_spec(snap) .* snap.masses
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



function get_x(snap::Snapshot)
    return get_x(snap.positions)
end

function get_x(A::Matrix{T}) where T<:Real
    return A[1, :]
end

function get_y(snap::Snapshot)
    return get_y(snap.positions)
end

function get_y(A::Matrix{T}) where T<:Real
    return A[2, :]
end

function get_z(snap::Snapshot)
    return get_z(snap.positions)
end

function get_z(A::Matrix{T}) where T<:Real
    return A[3, :]
end

function get_v_x(snap::Snapshot)
    return get_x(snap.velocities)
end

function get_v_y(snap::Snapshot)
    return get_y(snap.velocities)
end

function get_v_z(snap::Snapshot)
    return get_z(snap.velocities)
end

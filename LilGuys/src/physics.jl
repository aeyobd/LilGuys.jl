using LinearAlgebra: ×

"""
The magnitude of a matrix of 3-vectors of shape 3×N
"""
function calc_r(x::AbstractMatrix{T}) where T<:Real
    if size(x, 1) != 3
        error("Matrix must have shape 3×N")
    end
    return _calc_r(x)
end


function calc_r(a::AbstractArray{T}, b::AbstractArray{T}) where T<:Real
    return calc_r(a .- b)
end

function calc_r(x::AbstractVector{T}) where T<:Real
    if length(x) != 3
        error("Vector must have length 3")
    end
    return _calc_r(x)
end


function _calc_r(x::AbstractVector{T}) where T<:Real
    r = sqrt.(sum(x.^2, dims=1))
    if length(size(r)) == 1
        return r[1]
    else
        return dropdims(r, dims=1)
    end
end

function _calc_r(x::AbstractMatrix{T}) where T<:Real
    r = sqrt.(sum(x.^2, dims=1))
    return r[1, :]
end



function calc_r(snap::Snapshot)
    return calc_r(snap.positions)
end


function calc_v(snap::Snapshot)
    return calc_r(snap.velocities)
end

function calc_v(snap::Snapshot, v0::AbstractVector{T}) where T<:Real
    return calc_r(snap.velocities .- v0)
end

function calc_E_spec_kin(snap::Snapshot, v0=zeros(3))
    v = calc_v(snap, v0)
    return 0.5 .* v.^2
end


"""
Specific energy of each particle of a snapshot
"""
function calc_E_spec(snap::Snapshot)
    if snap.Φs == nothing || all(snap.Φs .== 0)
        Φ = calc_radial_discrete_Φ(snap)
    else
        Φ = snap.Φs
    end
    return calc_E_spec_kin(snap) .+ Φ
end


"""
    calc_ϵ(snap)

calculates binding energy of each particle in snapshot
"""
function calc_ϵ(snap::Snapshot)
    return -calc_E_spec(snap)
end


"""
Given potential and velocity, calculate specific energy
"""
function calc_E_spec(Φ::Real, v::Real)
    return 0.5v^2 .+ Φ
end


function calc_E_tot(snap::Snapshot, v_vec_0=zeros(3))
    if snap.Φs_ext == nothing
        Φs_ext = 0
    else
        Φs_ext = snap.Φs_ext
    end
    return sum(snap.masses .* (
               calc_E_spec_kin(snap, v_vec_0) 
               .+ 1/2 * snap.Φs 
               .+ Φs_ext)
              )
end


function calc_L_spec(snap::Snapshot)
    return calc_L_spec(snap.positions, snap.velocities)
end

function calc_L_spec(x::AbstractVector{T}, v::AbstractVector{T}) where T<:Real
    return x × v
end

function calc_L_spec(x::AbstractMatrix{T}, v::AbstractMatrix{T}) where T<:Real
    if size(x, 1) != 3 || size(v, 1) != 3
        error("Matrices must have shape 3×N")
    end
    if size(x, 2) != size(v, 2)
        error("Matrices must have same number of columns")
    end

    L = Matrix{F}(undef, 3, size(x, 2))

    for i in 1:size(x, 2)
        L[:, i] .= x[:, i] × v[:, i]
    end

    return L
end

function calc_L(snap::Snapshot)
    m = reshape(snap.masses, 1, :)
    return calc_L_spec(snap) .* m
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

function get_x(A::AbstractMatrix{T}) where T<:Real
    return A[1, :]
end

function get_y(snap::Snapshot)
    return get_y(snap.positions)
end

function get_y(A::AbstractMatrix{T}) where T<:Real
    return A[2, :]
end

function get_z(snap::Snapshot)
    return get_z(snap.positions)
end

function get_z(A::AbstractMatrix{T}) where T<:Real
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




"""
Returns a list of observations based on snapshot particles
"""
function to_sky(snap::Snapshot; invert_velocity::Bool=false, verbose::Bool=false)
    observations = Observation[]

    for i in 1:length(snap)
        if verbose
            print("converting $(i)/($(length(snap))\r")
        end
        pos = snap.positions[:, i]
        vel = snap.velocities[:, i]
        if invert_velocity
            vel *=-1
        end
        phase = PhasePoint(pos, vel)
        obs = to_sky(phase)
        push!(observations, obs)
    end
    return observations
end


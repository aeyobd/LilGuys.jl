using LinearAlgebra: ×



"""
    calc_r(x)

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



"""
    calc_r(snap, x_cen=zeros(3))

Calculates the radii of particles in a snapshot (from the center x_cen)
"""
function calc_r(snap::Snapshot, x_cen::AbstractVector{T}=snap.x_cen) where T<:Real
    return calc_r(snap.positions, x_cen)
end

"""The velocity of each particle in a snapshot"""
function calc_v(snap::Snapshot, v_cen::AbstractVector{T}=snap.v_cen) where T<:Real
    return calc_r(snap.velocities, v_cen)
end



get_x(snap::Snapshot) = get_x(snap.positions)
get_y(snap::Snapshot) = get_y(snap.positions)
get_z(snap::Snapshot) = get_z(snap.positions)

get_x(A::AbstractMatrix{T}) where T<:Real = A[1, :]
get_y(A::AbstractMatrix{T}) where T<:Real = A[2, :]
get_z(A::AbstractMatrix{T}) where T<:Real = A[3, :]

get_v_x(snap::Snapshot) = get_x(snap.velocities)
get_v_y(snap::Snapshot) = get_y(snap.velocities)
get_v_z(snap::Snapshot) = get_z(snap.velocities)




function calc_E_spec_kin(snap::Snapshot, v_cen=snap.v_cen)
    v = calc_v(snap, v_cen)
    return 0.5 .* v.^2
end


"""
    calc_E_spec(snap)

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


function calc_E_tot(snap::Snapshot, v_cen=snap.v_cen)
    if snap.Φs_ext == nothing
        Φs_ext = 0
    else
        Φs_ext = snap.Φs_ext
    end
    return sum(snap.masses .* (
               calc_E_spec_kin(snap, v_cen) 
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


"""
    calc_V_circ(snap, x_cen=zeros(3))

Returns a list of the sorted radii and circular velocity from a snapshot for the given centre.
"""
function calc_V_circ(snap::Snapshot, x_cen=snap.x_cen)
    r = calc_r(snap.positions .- x_cen)
    m = snap.masses[sortperm(r)]
    r = sort(r)
    M = cumsum(m)
    return r, calc_V_circ.(M, r)
end



"""
    get_bound(snap)

Returns the bound particles of a snapshot
"""
function get_bound(snap::Snapshot)
    return snap[E_spec(snap) .< 0]
end




"""
    to_sky(phase::PhasePoint)

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
        gc = Galactocentric(pos, vel)
        obs = transform(Observation, gc)
        push!(observations, obs)
    end
    return observations
end


import LinearAlgebra: ×


"""
    calc_r(x[, y])

The magnitude of a 3-vector or each vector in a matrix. Or, the distance between vecotrs x and y.
"""
function calc_r(x::AbstractMatrix{T}) where T<:Real
    if size(x, 1) != 3
        throw(DimensionMismatch("matrix must have 3 rows"))
    end
    return _calc_r(x)
end


function calc_r(a::AbstractArray{T}, b::AbstractArray{T}) where T<:Real
    return calc_r(a .- b)
end


function calc_r(x::AbstractVector{T}) where T<:Real
    if length(x) != 3
        throw(DimensionMismatch("Vector must have length 3"))
    end
    return _calc_r(x)
end


function _calc_r(x::AbstractVector{T}) where T<:Real
    r = sqrt.(sum(x.^2, dims=1))
    return r[1]
end


function _calc_r(x::AbstractMatrix{T}) where T<:Real
    r = sqrt.(sum(x.^2, dims=1))
    return r[1, :]
end



"""
    calc_r(snap[, x_cen])

Calculates the radii of particles in a snapshot (from the center x_cen).
Is stored in snapshot because very common calculation.
"""
function calc_r(snap::Snapshot, x_cen::AbstractVector{T}=snap.x_cen; recalculate=false) where T<:Real
    if snap._radii == nothing || recalculate
        snap._radii = calc_r(snap.positions .- x_cen)
    end
    return snap._radii
end



"""The velocity of each particle in a snapshot"""
function calc_v(snap::Snapshot, v_cen::AbstractVector{T}=snap.v_cen) where T<:Real
    return calc_r(snap.velocities, v_cen)
end



get_x(snap::Snapshot) = get_x(snap.positions)
get_y(snap::Snapshot) = get_y(snap.positions)
get_z(snap::Snapshot) = get_z(snap.positions)

get_x(A::AbstractMatrix{<:Real}) = A[1, :]
get_y(A::AbstractMatrix{<:Real}) = A[2, :]
get_z(A::AbstractMatrix{<:Real}) = A[3, :]

get_v_x(snap::Snapshot) = get_x(snap.velocities)
get_v_y(snap::Snapshot) = get_y(snap.velocities)
get_v_z(snap::Snapshot) = get_z(snap.velocities)




"""
    calc_K_spec(snap)

Kinetic energy of each particle of a snapshot
"""
function calc_K_spec(snap::Snapshot, v_cen=snap.v_cen)
    v2 = sum((snap.velocities .- v_cen) .^ 2, dims=1)
    return 1/2 * dropdims(v2, dims=1)
end


"""
    calc_K_tot(snap)

Total kinetic energy of a snapshot
"""
function calc_K_tot(snap::Snapshot)
    return sum(snap.masses .* calc_K_spec(snap))
end



@doc raw"""
    calc_W_tot(snap)

Total potential energy of a snapshot
```math
    W = -\frac{1}{2} \sum m_i Φ_i
```
"""
function calc_W_tot(snapshot::Snapshot)
    return -1/2 * sum(snapshot.masses .* snapshot.Φs)
end



"""
    calc_E_spec(snap)

Specific energy of each particle of a snapshot
"""
function calc_E_spec(snap::Snapshot)
    if snap.Φs == nothing 
        @warn "Snapshot does not contain Φ using radial Φ approximation."
        Φ = calc_radial_discrete_Φ(snap)
    else
        Φ = snap.Φs
    end
    return calc_K_spec(snap) .+ Φ
end


@doc raw"""
    calc_ϵ(snap)

calculates binding energy of each particle in snapshot.
```math
    ϵ = -E_{\text{spec}} = -\frac{1}{2}v^2 - Φ
```
"""
function calc_ϵ(snap::Snapshot)
    return -calc_E_spec(snap)
end



"""
    calc_E_spec(Φ, v)
Given potential and velocity, calculate specific energy.
"""
function calc_E_spec(Φ::Real, v::Real)
    return 0.5v^2 .+ Φ
end


"""
    calc_E_tot(snap)

Total energy of a snapshot
"""
function calc_E_tot(snap::Snapshot, v_cen=snap.v_cen)
    if snap.Φs_ext == nothing
        Φs_ext = 0
    else
        Φs_ext = snap.Φs_ext
    end
    return sum(snap.masses .* (
               calc_K_spec(snap, v_cen) 
               .+ 1/2 * snap.Φs 
               .+ Φs_ext)
              )
end



"""
    calc_L_spec(x, v)

Calculates the angular momentum of a particle with position x and velocity v
May pass a snapshot, or two 3-vecotrs, or two 3xN matrices for x and v.
"""
function calc_L_spec(x::AbstractVector{T}, v::AbstractVector{T}) where T<:Real
    return x × v
end


function calc_L_spec(snap::Snapshot)
    return calc_L_spec(snap.positions, snap.velocities)
end


function calc_L_spec(x::AbstractMatrix{T}, v::AbstractMatrix{T}) where T<:Real
    if size(x, 1) != 3 || size(v, 1) != 3
        throw(DimensionMismatch("Matrices must have 3 rows"))
    end
    if size(x, 2) != size(v, 2)
        throw(DimensionMismatch("Matrices must have the same number of columns"))
    end

    L = Matrix{F}(undef, 3, size(x, 2))

    for i in 1:size(x, 2)
        L[:, i] .= x[:, i] × v[:, i]
    end

    return L
end




"""
    calc_L_tot(snap)

Calculates the total angular momentum of a snapshot
"""
function calc_L_tot(snap::Snapshot)
    L = zeros(3)
    for i in 1:length(snap)
        L += snap.masses[i] .* calc_L_spec(snap.positions[:, i], snap.velocities[:, i])
    end

    return L
end


"""
    calc_v_circ(r, M)


The circular velocity at radius r from the center of a mass M.
"""
function calc_v_circ(r::Real, M::Real)
    if M < 0 || r < 0
        throw(DomainError("M and r must be positive"))
    elseif r == 0
        return 0
    end
    return  sqrt(G*M/r)
end


"""
    get_bound(snap)

Returns a filter for particles that are bound to the snapshot.
"""
function get_bound(snap::Snapshot)
    return calc_E_spec(snap) .< 0
end


"""
    get_bound_recursive_1D(snap)

Returns a filter for particles that are bound to the snapshot,
recursively removing unbound particles and recomputing the potential
"""
function get_bound_recursive_1D(snap::Snapshot; maxiter=300)
    r = calc_r(snap)
    m = snap.masses
    v = calc_v(snap)

    ϕ = calc_radial_discrete_Φ(r, m)
    ϵ = @. -1/2 * v^2 - ϕ
    filt = ϵ .> 0
    dN = sum(.!filt)

    for i in 1:maxiter
        if sum(filt) == 0
            break
        end

        ϕ = calc_radial_discrete_Φ(r[filt], m[filt])
        ϵ = @. -1/2 * v[filt]^2 - ϕ
        filt_2 = ϵ .> 0

        idx_dropped = eachindex(filt)[filt][.!filt_2]
        filt[idx_dropped] .= false

        dN = length(idx_dropped)

        if dN == 0
            break
        end
        if i == maxiter
            @warn "Maximum iterations reached"
        end
    end

    return filt
end


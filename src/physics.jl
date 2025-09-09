import LinearAlgebra: ×


"""
    radii(x[, y])

The magnitude of a 3-vector or each vector in a matrix. Or, the distance between vecotrs x and y.
"""
function radii(x::AbstractMatrix{T}) where T<:Real
    return _radii(x)
end


function radii(a::AbstractArray{T}, b::AbstractArray{T}) where T<:Real
    return radii(a .- b)
end


function radii(x::AbstractVector{T}) where T<:Real
    return _radii(x)
end


function _radii(x::AbstractVector{T}) where T<:Real
    r = sqrt.(sum(x.^2, dims=1))
    return r[1]
end


function _radii(x::AbstractMatrix{T}) where T<:Real
    r = sqrt.(sum(x.^2, dims=1))
    return r[1, :]
end



"""
    radii(snap[, x_cen])

Calculates the radii of particles in a snapshot (from the center x_cen).
Is stored in snapshot because very common calculation.
"""
function radii(snap::Snapshot, x_cen::AbstractVector{T}=snap.x_cen; recalculate=false) where T<:Real
    if snap._radii == nothing || recalculate
        snap._radii = radii(snap.positions .- x_cen)
    end
    return snap._radii
end



"""The velocity of each particle in a snapshot"""
function speeds(snap::Snapshot, v_cen::AbstractVector{T}=snap.v_cen) where T<:Real
    return radii(snap.velocities, v_cen)
end



x_position(snap::Snapshot) = x_position(snap.positions)
y_position(snap::Snapshot) = y_position(snap.positions)
z_position(snap::Snapshot) = z_position(snap.positions)

x_position(A::AbstractMatrix{<:Real}) = A[1, :]
y_position(A::AbstractMatrix{<:Real}) = A[2, :]
z_position(A::AbstractMatrix{<:Real}) = A[3, :]

x_velocity(snap::Snapshot) = x_position(snap.velocities)
y_velocity(snap::Snapshot) = y_position(snap.velocities)
z_velocity(snap::Snapshot) = z_position(snap.velocities)




"""
    kinetic_spec(snap)

Kinetic energy of each particle of a snapshot
"""
function kinetic_spec(snap::Snapshot, v_cen=snap.v_cen)
    v2 = sum((snap.velocities .- v_cen) .^ 2, dims=1)
    return 1/2 * dropdims(v2, dims=1)
end


"""
    kinetic_energy(snap)

Total kinetic energy of a snapshot
"""
function kinetic_energy(snap::Snapshot)
    return sum(snap.masses .* kinetic_spec(snap))
end



@doc raw"""
    potential_energy(snap)

Total potential energy of a snapshot
```math
    W = \frac{1}{2} \sum m_i Φ_i
```
"""
function potential_energy(snapshot::Snapshot)
    return 1/2 * sum(snapshot.masses .* snapshot.potential)
end



"""
    specific_energy(snap)

Specific energy of each particle of a snapshot (ϵ) defined to be positive for bound particles.
"""
function specific_energy(snap::Snapshot)
    if snap.potential == nothing 
        @warn "Snapshot does not contain Φ using radial Φ approximation."
        Φ = potential_spherical_discrete(snap)
    else
        Φ = snap.potential
    end
    return -(kinetic_spec(snap) .+ Φ)
end




"""
    specific_energy(Φ, v)
Given potential and velocity, calculate specific energy.
"""
function specific_energy(Φ::Real, v::Real)
    return -(1/2 * v^2 .+ Φ)
end


"""
    total_energy(snap)

Total energy of a snapshot
"""
function total_energy(snap::Snapshot, v_cen=snap.v_cen)
    if snap.potential_ext == nothing
        potential_ext = 0
    else
        Φs_ext = snap.potential_ext
    end
    return sum(snap.masses .* (
               kinetic_spec(snap, v_cen) 
               .+ 1/2 * snap.potential 
               .+ potential_ext)
              )
end



"""
    angular_momenta(x, v)

Calculate the angular momentum of a particle with position x and velocity v
May pass a snapshot, or two 3-vecotrs, or two 3xN matrices for x and v.
"""
function angular_momenta(x::AbstractVector{T}, v::AbstractVector{T}) where T<:Real
    return x × v
end


function angular_momenta(snap::Snapshot)
    return angular_momenta(snap.positions, snap.velocities)
end


function angular_momenta(x::AbstractMatrix{T}, v::AbstractMatrix{T}) where T<:Real
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
    angular_momentum(snap)

Calculates the total angular momentum of a snapshot
"""
function angular_momentum(snap::Snapshot)
    L = zeros(3)
    for i in 1:length(snap)
        L += snap.masses[i] .* angular_momenta(snap.positions[:, i], snap.velocities[:, i])
    end

    return L
end


"""
    v_circ(r, M)

The circular velocity at radius r from the center of a mass M.
"""
function v_circ(r::Real, M::Real)
    if M < 0 || r < 0
        throw(DomainError("M and r must be positive"))
    elseif r == 0
        return 0
    end
    return  sqrt(G*M/r)
end


"""
    bound_particles(snap)

Return a filter for particles that are bound to the snapshot.
"""
function bound_particles(snap::Snapshot; method=:simple, kwargs...)
    if method == :simple
        return specific_energy(snap) .>= 0
    elseif method == :nbody
        ϵ = -potential_nbody(snap) .- 1/2 .* speeds(snap).^2
        return ϵ .>= 0
    elseif method == :recursive_1D
        return bound_particles_recursive_1D(snap; kwargs...)
    elseif method == :recursive_3D
        return bound_particles_recursive_3D(snap; kwargs...)
    else
        throw(ArgumentError("Method $method not known"))
    end
end


"""
    bound_particles_recursive_1D(snap)

Find particles which are bound assuming spherical symmetry and
recursively removing unbound particles and updating the potential.
"""
function bound_particles_recursive_1D(snap::Snapshot; maxiter=300)
    r = radii(snap)
    m = snap.masses
    v = speeds(snap)

    ϕ = potential_spherical_discrete(r, m)
    ϵ = @. -1/2 * v^2 - ϕ
    filt = ϵ .> 0
    dN = sum(.!filt)

    for i in 1:maxiter
        @debug "bound particles iteration $(i-1), dutting $dN particles"
        if sum(filt) == 0
            break
        end
        ϕ = potential_spherical_discrete(r[filt], m[filt])
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


"""
    bound_particles_recursive_3D(snap)

Find particles which are bound recursively removing unbound particles and updating the potential.
"""
function bound_particles_recursive_3D(snap::Snapshot; maxiter=300)
    v = speeds(snap)

    ϕ = potential_nbody(snap)
    ϵ = @. -1/2 * v^2 - ϕ
    filt = ϵ .> 0
    dN = sum(.!filt)

    for i in 1:maxiter
        if sum(filt) == 0
            break
        end
        snap_new = snap[filt]

        ϕ = potential_nbody(snap_new)
        ϵ = @. -1/2 * v[filt]^2 - ϕ
        filt_2 = ϵ .> 0

        idx_dropped = eachindex(filt)[filt][.!filt_2]
        filt[idx_dropped] .= false

        dN = length(idx_dropped)

        @debug "bound particles iteration $i, dutting $dN particles"
        if dN == 0
            break
        end
        if i == maxiter
            @warn "Maximum iterations reached"
        end
    end

    return filt
end


""" Abstract type representing a force modifying particles in an orbit"""
abstract type Potential end

function potential(p::Potential, pos)
    error("calc_Φ not implemented for this potential")
end

function acceleration(p::Potential, pos, vel=nothing, t=nothing)
    error("calc_Φ not implemented for this potential")
end


struct CompositePotential <: Potential
    potentials::Vector{Potential}
end

function potential(pot::CompositePotential, pos::AbstractVector{<:Real})
    return sum(potential(p, pos) for p in pot.potentials)
end
function acceleration(pot::CompositePotential, pos, vel, t)
    return sum(acceleration(p, pos, vel, t) for p in pot.potentials)
end

Base.@kwdef struct ChandrashakarDynamicalFriction <: Potential
    r_s::Real
    σv::Function
    ρ::Function
    M::Real
    Λ::Union{Real, Nothing} = nothing
    ϵ::Real = ϵ_default(r_s)
end


function potential(p::ChandrashakarDynamicalFriction, pos)
    return 0.0
end


function acceleration(p::ChandrashakarDynamicalFriction, pos, vel, t=nothing)
    return a_dyn_friction(pos, vel; r_s=p.r_s, σv=p.σv, ρ=p.ρ, M=p.M, Λ=p.Λ, ϵ=p.ϵ)
end


@doc raw"""
    a_dyn_friction(pos, vel; σv, ρ, M, [Λ, r_s])

Compute the Chandrashakar dynamical friction.

``
\frac{d{\bf v}}{dt} = -\frac{4\pi\,G^2\,M\,\rho\,\ln\Lambda}{v^2} \left({\rm erf}(X) - \frac{2X}{\sqrt\pi} \exp(-X^2)\right) \frac{{\bf v}}{v}
``

Assumes the host is at the origin, σv and ρ are functions of a 3-vector position giving velocity dispersion and density respectively, M is the satellite mass, and Λ is the coloumb integral.
If the coloumb integral is not given, then the satellite scale radius should be provided (`r_s`) and Λ is approximated as 

``
\Lambda = \begin{cases}
    r / 0.45r_s & r_s < 8 \\ 
    r / (2.2r_s - 14) & {\rm otherwise}
\end{cases}
``

where $r$ is the current radius of the satellite.
"""
function a_dyn_friction(pos, vel; σv, ρ, M, Λ=nothing, r_s=nothing, ϵ=ϵ_default(r_s))
    v = radii(vel)
    
    X = v / (√2 * σv(pos))
    fX = erf(X) - 2X/√π * exp(-X^2)

    r = LilGuys.radii(pos)
    if isnothing(Λ)
        Λ = Λ_default(r, ϵ)
    end

    return -4π*G^2 * M * ρ(pos) * max(log(Λ), 0) * fX * vel ./ v^3
end


function ϵ_default(r_s)
    if r_s < 8
        ϵ = 0.45r_s
    else
        ϵ = 2.2r_s - 14
    end
    return ϵ
end


function Λ_default(r, ϵ)
    Λ = r / ϵ
    return max(Λ, 1) # should not be less than 1
end



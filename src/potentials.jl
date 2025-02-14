
""" Abstract type representing a force modifying particles in an orbit"""
abstract type Potential end

function update_potential!(p::Potential, t::Real)
    return
end

function calc_Φ(p::Potential, pos::Vector{F})
    error("calc_Φ not implemented for this potential")
end

function calc_a(p::Potential, pos::Vector{F})
end


struct CompositePotential <: Potential
    potentials::Vector{Potential}
end


function update_potential!(p::CompositePotential, t::Real)
    for pot in p.potentials
        update_potential!(pot, t)
    end
end


struct ChandrashakarDynamicalFriction <: Potential
    r_s::Real
    σv::Function
    ρ::Function
    M::Real
    Λ::Real
end

function calc_Φ(p::ChandrashakarDynamicalFriction, pos::Vector{F})
    return 0.0
end

function calc_a(p::ChandrashakarDynamicalFriction, pos::Vector{F})
    return a_dyn_friction(pos, vel; r_s=p.r_s, σv=p.σv, ρ=p.ρ, M=p.M, Λ=p.Λ)
end


@doc raw"""
    a_dyn_friction(pos, vel; r_s, σv, ρ, M, Λ)

Computes the dynamical friction based on 

``
\frac{d{\bf v}}{dt} = -\frac{4\pi\,G^2\,M\,\rho\,\ln\Lambda}{v^2} \left({\rm erf}(X) - \frac{2X}{\sqrt\pi} \exp(-X^2)\right) \frac{{\bf v}}{v}
``
"""
function a_dyn_friction(pos, vel; r_s=nothing, σv, ρ, M, Λ=nothing  )
    G = 1
    v = calc_r(vel)
    r = calc_r(pos)
    
    if Λ === nothing
        if r_s < 8
            ϵ = 0.45r_s
        else
            ϵ = 2.2r_s - 14
        end

        Λ = r / ϵ
    end

    X = v / (√2 * σv(r))
    fX = erf(X) - 2X/√π * exp(-X^2)

    return -4π*G^2 * M * ρ(pos) * max(log(Λ), 0) * fX * vel ./ v^3
end


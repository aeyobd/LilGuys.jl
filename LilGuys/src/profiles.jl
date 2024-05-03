import QuadGK: quadgk 
"""A series of methods describing common mass profiles"""
import Base: @kwdef

abstract type AbstractProfile


struct ABC
    α::Float64
    β::Float64
    γ::Float64
    M::Float64 = 1
    r_s::Float64 = 1
end


struct Plummer
    a::Float64
    M::Float64 = 1
    r_s::Float64 = 1
end


struct Sersic
end



struct Exp2D
    M::Float64 = 1
    R_s::Float64 = 1
end


struct Exp3D
    M::Float64 = 1
    r_s::Float64 = 1
end


struct NFW
    M_s::Float64 = 1
    r_s::Float64 = 1
end



function calc_ρ(profile::ABC_Profile, r::Float64)
    α, β, γ = profile.α, profile.β, profile.γ
    return r^(-γ) * (1 + r^a) ^ ((γ - β)/α)
end


function calc_ρ1(profile::ABC_Profile, r::Float64)
    α, β, γ = profile.α, profile.β, profile.γ
    return r^(-γ) * (1 + r^a) ^ ((γ - β)/α)
end



ρ_s_int(r, r_s, n) = -r_s^n * exp(-(r/r_s)^n) * (r/r_s)^n * (n+1)

function calc_ρ(profile::Exp3D, r::Float64)
    M, r_s = profile.M, profile.r_s
    return M / (8π * r_s^3) * exp(-r/r_s)
end


function M(profile::Exp3D, r::Float64)
    M, r_s = profile.M, profile.r_s
    return 
end




function Σ(profile::Exp2D, r::Float64)
    M, R_s = profile.M, profile.R_s
    return M / (2π * R_s^2) * exp(-r/R_s)
end




function V_circ(profile, r)
    return sqrt(G * M(profile, r) / r)
end

function M(profile, r::Float64)
    return quadgk(r -> 4π * r^2 * calc_ρ(profile, r), 0, r)[1]
end

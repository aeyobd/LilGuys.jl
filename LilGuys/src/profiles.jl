import QuadGK: quadgk 
"""A series of methods describing common mass profiles"""

abstract type AbstractProfile

# check these!
ρ_s(r, r_s, n) = exp(-(r/r_s)^n)
ρ_s_int(r, r_s, n) = -r_s^n * exp(-(r/r_s)^n) * (r/r_s)^n * (n+1)

struct ABC
    α::Float64
    β::Float64
    γ::Float64
    M::Float64 = 1
end


struct Plummer
    a::Float64
    M::Float64 = 1
end


struct Sersic
end


function calc_ρ(profile::ABC_Profile, r::Float64)
    α, β, γ = profile.α, profile.β, profile.γ
    return r^(-γ) * (1 + r^a) ^ ((γ - β)/α)
end

function calc_ρ1(profile::ABC_Profile, r::Float64)
    α, β, γ = profile.α, profile.β, profile.γ
    return r^(-γ) * (1 + r^a) ^ ((γ - β)/α)
end

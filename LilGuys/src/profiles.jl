import QuadGK: quadgk 


struct ABC_Profile
    α::Float64
    β::Float64
    γ::Float64
    M::Float64 = 1
end


struct PlummerProfile
    a::Float64
    M::Float64 = 1
end


function calc_ρ(profile::ABC_Profile, r::Float64)
    α, β, γ = profile.α, profile.β, profile.γ
    return r^(-γ) * (1 + r^a) ^ ((γ - β)/α)
end

function calc_ρ1(profile::ABC_Profile, r::Float64)
    α, β, γ = profile.α, profile.β, profile.γ
    return r^(-γ) * (1 + r^a) ^ ((γ - β)/α)
end

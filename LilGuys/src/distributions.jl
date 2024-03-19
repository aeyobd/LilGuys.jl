import QuadGK: quadgk


"""
    calc_fϵ(ρ, ψ; n_r=100, n_E=100, r_min=1e-2, r_max=100)

given functions ρ(r) and ψ(r), computes the binding energy 
distribution function f(ϵ) using the Eddington inversion formula.
"""
function calc_fϵ(ρ::AbstractVector, ψ::AbstractVector, r::AbstractVector)
    if !issorted(r) | !issorted(-ψ)
        throw(ArgumentError("arrays must be sorted"))
    end

    ρ1 = gradient(ρ, r)
    ρ2 = gradient(ρ1, r)
    ψ1 = gradient(ψ, r)
    ψ2 = gradient(ψ1, r)

    d2ρ_dψ2 = @. ψ1^-2 * ρ2 - ψ1^-3 * ρ1 * ψ2
    d2_interp = lerp(reverse(ψ), reverse(d2ρ_dψ2))

    return calc_fϵ_from_∂(d2_interp)
end


"""
    calc_fϵ_from_∂(ψ, d2ν_dψ2)

given (interpolated) function d2ν_dψ2(ψ), compute f(ϵ) from the Eddington inversion formula.
"""
function calc_fϵ_from_∂(d2ν_dψ2)
    f_integrand(ϵ, ψ) = 1/√8π^2 * d2ν_dψ2(ψ) /√(ϵ - ψ)
    return ϵ -> quadgk(ψ->f_integrand(ϵ, ψ), 0, ϵ)[1]
end


"""
    warn_missed_regions(ρ, r_bins)

Utility function to print out the fraction of the mass that is outside binned range. Works for a callable ρ
"""
function warn_missed_regions(ρ, r_bins)
    M_integrand = r -> 4π * r^2  * ρ(r)
    M(r) = quadgk(M_integrand, 0, r)[1]
    M_tot = M(Inf)
    
    # print out excluded particles due to binning
    fin = M(r_bins[1])
    fout = M_tot - M(r_bins[end])
    println("fraction ρ inside bins $fin, fraction outside $fout")
end

#!/usr/bin/env julia
using ArgParse
import QuadGK: quadgk
import LilGuys as lguys

function get_args()
    # Parse command line arguments
    s = ArgParseSettings(description = "generates an exopential n-body equilibrium profile. Uses that 1R = 1 half light radius in 2D, and mass integrates to 1. ")

    @add_arg_table s begin
        "output"
            help="Output file"
            default="model.hdf5"
        "-m", "--mass"
            help="total stellar mass in 1e10 solar mass units"
            default=1e-6
            arg_type=Float64
        "-r", "--radius"
            help="the scale radius"
            arg_type=Float64

        "--nr"
            help="Number of radial points"
            default=10000
        "--nE"
            help="Number of energy points"
            default=1000
        "--rmin"
            help="Minimum radius"
            default=1e-2
        "--rmax"
            help="Maximum radius"
            default=1e2
    end


    return parse_args(s)
end


function ρ_exp(r)
    return 1/π * exp(-2r)
end


function kroupa_imf(m::Real)
    α0 = 0.3
    α1 = 1.3
    α2 = 2.3
    α3 = 2.7

    mmin = 0.01
    m01 = 0.08
    m12 = 0.50
    m23 = 1

    k2 = 1 / 7.7279
    k1 = k2 * m12^(α1-α2)
    k0 = k1 * m01^(α0-α1)
    k3 = k2 * m23^(α3-α2)

    if mmin <=m < m01
        return k0 * m^-α0
    elseif m01 <= m < m12
        return k1 * m^-α1
    elseif m12 <= m < m23
        return k2 * m^-α2
    elseif m23 <=m
        return k3 * m^-α3
    else
        error("only accepts real input, got $m")
    end

end


function calc_m_tot(imf, m_min, m_max)
    return quadgk(m->imf(m) * m, m_min, m_max)
end

function calc_n_tot(imf, m_min, m_max)
    return quadgk(imf, m_min, m_max)
end

function calc_m_mean(imf=kroupa_imf, m_min=0.08, m_max=100)
    return calc_m_tot(imf, m_min, m_max)[1] / calc_n_tot(imf, m_min, m_max)[1]
end

function calc_inv_imf(;imf=kroupa_imf, m_res=100_000, m_min=0.08, m_max=100)
    ms = 10 .^ LinRange(log10(m_min), log10(m_max), m_res)
    dm = lguys.gradient(ms)
    imfs = imf.(ms)
    probs = imfs .* dm

    cdf = cumsum(probs)
    cdf ./= cdf[end]
    return lguys.lerp(cdf, ms)
end


function sample_ms(;imf=kroupa_imf, m_min=0.08, m_max=100, m_tot=1e-6, m_res=10_000)
    m_from_p = calc_inv_imf(;m_min=m_min, m_max=m_max, m_res=m_res, imf=imf)

    m_mean = calc_m_mean(imf, m_min, m_max)
    n_samples = lguys.M0 * m_tot / m_mean
    println(n_samples)

    ps = lguys.randu(0, 1, round(Int64, n_samples))
    masses = m_from_p.(ps)

    return masses
end




function main(ρ_in=ρ_exp)
    G = 1
    args = get_args()
    NE = round(Int, args["nE"])
    r_bins = 10 .^ range(log10(args["rmin"]), log10(args["rmax"]), length=args["nr"])
    println("sampling masses")
    masses = sample_ms(m_tot = args["mass"])
    masses ./= lguys.M0
    println("mass total = $(sum(masses))")
    N = length(masses)
    println("n particles = $N")

    Mtot = 4π * quadgk(x->x^2 * ρ_in(x), 0, Inf)[1]
    
    ρ(r) = ρ_in(r) / Mtot
    M(r) = 4π* quadgk(x->x^2 * ρ(x), 0, r)[1]

    radii = sample_radii(M, r_bins, N)

    Φ(r) = -G*M(r)/r  - 4π*G * quadgk(x->x * ρ(x), r, Inf)[1]
    ψs = -Φ.(r_bins)
    ν = ρ.(r_bins)
    Emin = min(ψs[end], ψs[end] / NE)
    E = LinRange(Emin, ψs[1], NE)

    ψ_i = lguys.lerp(r_bins, ψs)
    ϵs = sample_ϵs(ν, ψ_i, r_bins, E, radii)
    vs = @. sqrt(2*(ψ_i(radii) - ϵs))

    x_hat = lguys.rand_unit(N)
    v_hat = lguys.rand_unit(N)
    positions = x_hat .* radii'
    positions *= args["radius"]
    velocities = v_hat .* vs'

    m_scale = sum(masses) 
    v_scale = sqrt(lguys.G * m_scale / args["radius"])
    velocities .*= v_scale

    snap = lguys.Snapshot(positions, velocities, masses)


    println("Saving snapshot to $(args["output"])")
    lguys.save(args["output"], snap)

end


function sample_ϵs(ν, ψ, r_bins, E, radii)
    P = P_ϵ_r(ν, ψ, r_bins, E)

    N = length(radii)
    ϵs = Vector{Float64}(undef, N)
    n_rej = 0
    for i in 1:N
        ϵs[i], n = sample_ϵ(radii[i], ψ.(radii[i]), P)
        n_rej += n
    end

    rate = N / (N + n_rej)
    println("acceptance rate = $rate")
    return ϵs
end


function sample_ϵ(r, ψ, P)
    n_rej = 0
    ϵ = ψ .* rand()

    while true
        ϵ = ψ .* rand()
        p = P(ϵ, r)
        n_rej += 1

        if rand() < p
            break
        end
    end

    return ϵ, n_rej
end


function P_ϵ_r(ν, ψ_i , r_bins, E)
    ψs = ψ_i.(r_bins)
    ψ_i = lguys.lerp(r_bins, ψs)
    f_i = f_interp(ν, ψs, r_bins, E)
    r_i = lguys.lerp(reverse(ψs), reverse(r_bins))

    dP_dr(ϵ, r) = sqrt(2*(ψ_i(r) - ϵ)) * r^2 
    L_P(ϵ, r) = f_i(ϵ) * dP_dr(ϵ, r) # likelihood function

    dP(ϵ) = quadgk(r->dP_dr(ϵ, r), 0, r_i(ϵ))[1]

    L_P_max = [maximum(L_P.(E[E .<= ψ_i(x)], x)) for x in r_bins]
    L_P_max .*= 1.1
    L_P_max_i = lguys.lerp(r_bins, L_P_max)

    return (ϵ, r) -> L_P(ϵ, r) / L_P_max_i(r)
end



function f_interp(ν, ψs, r_bins, E)
    f = lguys.calc_fϵ(ν, ψs, r_bins)
    
    fE = f.(E)

    if any(fE .< 0)
        println("$(sum(fE .< 0)) f are  < 0")
        fE[fE .< 0] .= 0
    end
    f_i = lguys.lerp(E, fE)
    return f_i
end


function sample_radii(M::Function, r_bins, N)
    Mtot = M(Inf)
    Min = M(r_bins[1])
    Mout = Mtot - M(r_bins[end])
    println("missing Min = $Min, Mout = $Mout")
    r_m_i = lguys.lerp(M.(r_bins), r_bins)
    m = lguys.randu(Min, Mtot-Mout, N)
    r = r_m_i.(m)
    return r
end




if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

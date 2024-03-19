using ArgParse
import QuadGK: quadgk
import LilGuys as lguys

function get_args()
    # Parse command line arguments
    s = ArgParseSettings()

    @add_arg_table s begin
        "output"
            help="Output file"
            default="model.hdf5"
        "-n"
            help="number of particles"
            default=10000
        "--nr"
            help="Number of radial points"
            default=100
        "--nE"
            help="Number of energy points"
            default=100
        "--rmin"
            help="Minimum radius"
            default=1e-2
        "--rmax"
            help="Maximum radius"
            default=1e2
    end


    return parse_args(s)
end

function main()
    args = get_args()

    ρ(r) = r^-γ * (1 + r^α)^((γ-β)/α)
    α = 1
    β = 4
    γ = 1
    G = 1

    M(r) = 4π * quadgk(x->x^2 * ρ(x), 0, r)[1]
    Φ(r) = -G*M(r)/r  - 4π*G * quadgk(x->x * ρ(x), r, Inf)[1]

    r = 10 .^ range(log10(args["rmin"]), log10(args["rmax"]), length=args["nr"])

    Mtot = M(Inf)
    ψ = -Φ.(r) ./ Mtot
    ν = ρ.(r) ./ Mtot

    f = lguys.calc_fϵ(ν, ψ, r)
    E = LinRange(ψ[end], ψ[1], args["nE"])
    fE = f.(E)

    if any(fE .< 0)
        println("$(sum(f .< 0)) f are  < 0")
        fE[fE .< 0] .= 0
    end

    Min = M(r[1])
    Mout = Mtot - M(r[end])
    println("Mtot = $Mtot, Min = $Min, Mout = $Mout")

    ψ_i = lguys.lerp(r, ψ)
    ν_i = lguys.lerp(r, ν)
    f_i = lguys.lerp(r, f.(E))
    r_i = lguys.lerp(reverse(ψ), reverse(r))
    r_m_i = lguys.lerp(M.(r), r)

    dP_dr(ϵ, r) = sqrt(2*(ψ_i(r) - ϵ)) * r^2 
    L_P(ϵ, r) = f_i(ϵ) * dP_dr(ϵ, r)

    dP(ϵ) = quadgk(r->dP_dr(ϵ, r), 0, r_i(ϵ))[1]

    L_P_max = [maximum(L_P.(E[E .<= ψ_i(x)], x)) for x in r]
    L_P_max .*= 1.1

    L_P_max_i = lguys.lerp(r, L_P_max)

    N = round(Int64, args["n"])
    radii = zeros(N)
    vels = zeros(N)

    n = 0
    k = 1
    n_rej = 0
    while n < N
        ms = rand(k) * (Mtot - Mout - Min) .+ Min
        rs = r_m_i.(ms)
        ψs = ψ_i.(rs)

        L_P_maxs = L_P_max_i.(rs)
        ϵs = ψs .* rand(k)
        Ls = L_P.(ϵs, rs)
        Ys = rand(k) .* L_P_maxs
        while Ys > Ls
            ϵs = ψs .* rand(k)
            Ls = L_P.(ϵs, rs)
            Ys = rand(k) .* L_P_maxs
            n_rej += 1
        end
        n += k
        radii[n:n+k-1] .= rs
        v = sqrt.(2 * (ψs .- ϵs))
        vels[n:n+k-1] .= v
    end

    x_hat = lguys.rand_unit(N)
    v_hat = lguys.rand_unit(N)
    pos = x_hat .* reshape(radii, (1, N))
    vel = v_hat .* reshape(vels, (1, N))
    masses = ones(N) / N

    snap = lguys.Snapshot(pos, vel, masses)

    lguys.save(args["output"], snap)

end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

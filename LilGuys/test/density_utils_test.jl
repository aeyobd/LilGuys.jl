
@testset "calc_Σ" begin
    log_r = [-Inf, 0]
    mass = [1]
    Σ = LilGuys.calc_Σ(log_r, mass)
    @test Σ ≈ [1/π]

    log_r = [-1, 0, 1, 1.5]
    mass = [1, 2, 0]
    Σ = LilGuys.calc_Σ(log_r, mass)
    @test Σ ≈ [1 / (0.99π), 2 / (99π), 0]
end


function sample_density_profile(f, n)
    log_r = LinRange(-5, 5, 1000)
    r = exp10.(log_r)
    Σ = f.(r)
    M = cumsum(Σ .* π .* r .* LilGuys.gradient(r))
    M ./= M[end]
    l = LilGuys.lerp([0; M], [0; r])

    probs = rand(n)
    return l.(probs)
end

@testset "integration with exp profile" begin
    # TODO: normalize this; incorperate predicted uncertainties
    f(r) = exp(-r)
    N = 10_000
    r = sample_density_profile(f, N)

    mass = 0.5 .+ 0.5rand(N)
    M = sum(mass)

    obs = LilGuys.calc_properties(r, normalization=:none, weights=mass)

    @test sum(obs.mass_in_annulus) ≈ sum(mass)
    @test obs.M_in[end] ≈ sum(mass)
    @test issorted(obs.log_r)
    @test sum(obs.counts) ≈ N
    @test obs.Sigma ≈ M * f.(10 .^ obs.log_r_bins[2:end])
    # TODO: Gamma, Gamma max, 

    # TODO: try other normalizations
end

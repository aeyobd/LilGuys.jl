test_profiles = [
    lguys.Exp2D(1, 1),
    lguys.Exp2D(1.3, 0.7124),
    lguys.Exp3D(1.3, 0.7124),
    lguys.Exp3D(1, 1),
    lguys.LogCusp2D(1, 1),
    lguys.KingProfile(M=0.52, R_s=1.123, R_t=2π),
]

dir = mktempdir()

import TOML


@testset "Exp2D Σ" begin
    M = 1
    r_s = 1
    profile = lguys.Exp2D(M, r_s)

    @test lguys.calc_Σ(profile, 0) ≈ 1 / 2π * M / r_s^2
    @test lguys.calc_Σ(profile, r_s) ≈ exp(-1) / 2π * M / r_s^2
    @test lguys.calc_Σ(profile, 2r_s) ≈ exp(-2) / 2π * M / r_s^2
    @test lguys.calc_Σ(profile, Inf) ≈ 0
end


@testset "Exp3D ρ" begin
    M = 1
    r_s = 1

    profile = lguys.Exp3D(M, r_s)
    
    ρ0 = 1 / 4π * M / r_s^3
    @test lguys.calc_ρ(profile, 0) ≈ ρ0 broken=true

end


@testset "log cusp 2D Σ" begin
    M = 1
    R_s = 1

    profile = lguys.LogCusp2D(M, R_s)

    @test lguys.calc_Σ(profile, 0) ≈ Inf
    @test lguys.calc_Σ(profile, R_s) ≈ 0 broken=true
end


@testset "King" begin
    M = 1
    r_s = 1
    r_t = 2
    profile = lguys.KingProfile(k=1, R_s=1, R_t=r_t)

    @test lguys.calc_Σ(profile, 0) ≈ (1 - (1+(r_t/r_s)^2)^(-1/2))^2

    @test lguys.calc_Σ(profile, r_t) ≈ 0
    @test lguys.calc_Σ(profile, 1.2*r_t) ≈ 0
    @test lguys.calc_Σ(profile, Inf) ≈ 0

    @test lguys.calc_M_2D(profile, r_t) ≈ profile.M rtol=1e-10 broken=true

end

@testset "total mass" begin
    for profile in test_profiles
        if profile isa lguys.KingProfile
            @test lguys.calc_M_2D(profile, 100.0) ≈ 0.52
        else
            @test lguys.calc_M(profile, 100.0) ≈ profile.M
        end
    end
end


@testset "total 2D mass" begin
    for profile in test_profiles
        integrand(r) = 2 * π * r * lguys.calc_Σ(profile, r)

        @test lguys.quadgk(integrand, 0, Inf)[1] ≈ profile.M
    end
end


@testset "3d to 2d density" begin

    for profile in test_profiles
        integrand(r, R) = 2 *r* lguys.calc_ρ(profile, r) / sqrt(r^2 - R^2)

        x = 10 .^ LinRange(-2, 0, 100)
        eps = 1e-6
        Σ(R) = lguys.quadgk(r -> integrand(r, R), R*(1+eps), Inf)[1]

        @test lguys.calc_Σ.(profile, x) ≈ Σ.(x) rtol=1e-3
    end
end


@testset "2d to 3d density" begin
    x = [0, 0.015, 0.23, 0.95, 1.254, 2.53, 3.0, 5., 100]
    for profile in test_profiles
        ρ1 = LilGuys.calc_ρ.(profile, x)
        ρ2 = LilGuys.calc_ρ_from_Σ.(profile, x)

        @test ρ1 ≈ ρ2 rtol=1e-5
    end
end


@testset "calc_R_h" begin
    for profile in test_profiles
        R_h = lguys.calc_R_h(profile)
        @test lguys.calc_M_2D(profile, R_h) ./ profile.M ≈ 0.5 rtol=1e-3
    end
end


@testset "calc_r_h" begin
    for profile in test_profiles
        r_h = lguys.calc_r_h(profile)
        @test lguys.calc_M(profile, r_h) ./ profile.M ≈ 0.5 rtol=1e-3
    end
end

@testset "load profile" begin

    @testset "basic tests" begin
        d = Dict(
            "Exp2D" => Dict("R_s" => 2.0, "M" => 1.5),
           )
        
        prof = lguys.load_profile(d)
        @test prof.R_s == 2.0
        @test prof.M == 1.5

        filename = joinpath(dir, "test_profile.toml")
        open(filename, "w") do io
            TOML.print(io, d)
        end

        prof = lguys.load_profile(filename)
        @test prof.R_s == 2.0
        @test prof.M == 1.5


        d = Dict(
                 "profile" => Dict("Exp2D" => Dict("R_s" => 0.1, "M" => 1.2))
           )
        
        prof = lguys.load_profile(d)
        @test prof.R_s == 0.1
        @test prof.M == 1.2

        open(filename, "w") do io
            TOML.print(io, d)
        end

        prof = lguys.load_profile(filename)
        @test prof.R_s == 0.1
        @test prof.M == 1.2
    end

    @testset "error handling" begin
        d = Dict{String, Any}()
        @test_throws ArgumentError lguys.load_profile(d)

        d = Dict("profile" => Dict{String, Any}())
        @test_throws ArgumentError lguys.load_profile(d)

        d = Dict("profile" => Dict("Exp2D" => Dict(), "Exp3d" => Dict()))
        @test_throws ArgumentError lguys.load_profile(d)


        d = Dict("profile" => Dict("Exp2D" => Dict("R_s" => 0.1, "M" => 1.2), "Exp3D" => Dict("R_s" => 0.1, "M" => 1.2)))
        @test_throws ArgumentError lguys.load_profile(d)

        d = Dict("profile" => Dict("GooblyGawk" => Dict("R_s" => 0.1, "M" => 1.2)))
        @test_throws KeyError lguys.load_profile(d)
    end

    @testset "retrieval" begin
        for profile in lguys.RECOGNIZED_PROFILES
            for _ in 1:1000
                kwargs = Dict()

                profile_class = getproperty(lguys, profile)

                for k in fieldnames(profile_class)
                    if rand() < 0.5
                        kwargs[string(k)] = rand()
                    end
                end

                if profile == :KingProfile
                    kwargs["R_t"] = 10rand()
                end

                if profile == :NFW
                    kwargs = Dict()
                    kwargs["M_s"] = 10^randn()
                    kwargs["r_s"] = 10rand()
                    kwargs["c"] = 10rand()
                end

                d = Dict("profile" => Dict(string(profile) => kwargs))
                prof = lguys.load_profile(d)

                prof_exp = profile_class(; lguys.dict_to_tuple(kwargs)...)

                for k in keys(d["profile"][string(profile)])
                    k = Symbol(k)
                    @test getfield(prof, k) == getfield(prof_exp, k)
                end

            end
        end
    end
end

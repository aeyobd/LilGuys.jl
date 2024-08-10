test_profiles = [
    lguys.Exp2D(1, 1),
    lguys.Exp2D(1.3, 0.7124),
    lguys.Exp3D(1.3, 0.7124),
    lguys.Exp3D(1, 1),
    lguys.LogCusp2D(1, 1),
]

import TOML

@testset "total mass" begin
    for profile in test_profiles
        @test lguys.calc_M(profile, 100.0) ≈ profile.M
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

        @test lguys.calc_Σ.([profile], x) ≈ Σ.(x) rtol=1e-3
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

        filename = "test_profile.json"
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

        filename = "test_profile.json"
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

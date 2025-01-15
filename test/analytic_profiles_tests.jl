dir = mktempdir()

import TOML

function test_R_h(profile, M=profile.M)
    R_h = lguys.calc_R_h(profile)
    @test lguys.calc_M_2D(profile, R_h) ./ M ≈ 1/2 rtol=1e-3
end

function test_r_h(profile, M=profile.M)
    r_h = lguys.calc_r_h(profile)
    @test lguys.calc_M(profile, r_h) ./ M ≈ 1/2 rtol=1e-3
end

function test_M_tot(profile, M=profile.M; r=100)
    @test lguys.calc_M(profile, r) ≈ M
end

function test_M_2D_tot(profile, M=profile.M, r=100)
    @test lguys.calc_M_2D(profile, r) ≈ M
end

function test_to_zero(profile; r=100000, atol=1e-8)
    @test lguys.calc_Σ(profile, r) ≈ 0 atol=atol
    @test lguys.calc_ρ(profile, r) ≈ 0 atol=atol
    @test lguys.calc_M(profile, 1/r) ≈ 0 atol=atol
end


@testset "Plummer" begin
    @testset "Σ" begin
        profile = lguys.Plummer(1, 1)
        @test lguys.calc_Σ(profile, 0) ≈ 1 / π
        @test lguys.calc_Σ(profile, 1) ≈ 1 / (4π)
        @test lguys.calc_Σ(profile, 2) ≈ 1 / (25π)

        @test lguys.calc_R_h(profile) ≈ 1
    end

    @testset "ρ" begin
        profile = lguys.Plummer(0.23, 0.51)
        x = 10 .^ LinRange(-1, 1, 10)
        ρ1 = lguys.calc_ρ.(profile, x)
        ρ2 = lguys.calc_ρ_from_Σ.(profile, x)
        @test ρ1 ≈ ρ2 rtol=1e-5
    end

    @testset "Σ_inv" begin
        profile = lguys.Plummer(0.23, 0.51)
        x = 10 .^ LinRange(-1, 1, 10)
        Σ1 = lguys.calc_Σ.(profile, x)
        Σ2 = lguys.calc_Σ_from_ρ.(profile, x)
        @test Σ1 ≈ Σ2 rtol=1e-5
    end

    @testset "M" begin
        profile = lguys.Plummer(0.989, 1.35)
        x = 10 .^ LinRange(-1, 0.5, 10)

        M1 = lguys.calc_M.(profile, x)
        M2 = lguys.calc_M_from_ρ.(profile, x)

        @test M1 ≈ M2 rtol=1e-5
    end


    @testset "M2D" begin
        profile = lguys.Plummer(0.989, 1.35)
        x = 10 .^ LinRange(-1, 0.5, 10)

        M1 = lguys.calc_M_2D.(profile, x)
        M2 = lguys.calc_M_2D_from_Σ.(profile, x)

        @test M1 ≈ M2 rtol=1e-5
    end

    @testset "consistency" begin
        profile = lguys.Plummer(0.987, 3.3)

        test_R_h(profile)
        test_r_h(profile)
        test_M_tot(profile, r=1e5)
        test_to_zero(profile)
    end

    @testset "scale" begin
        prof = lguys.Plummer(2.3, 0.6)
        r_scale = 1.21
        m_scale = π

        prof_scaled = lguys.scale(prof, r_scale, m_scale)

        x = 10 .^ LinRange(-1, 1, 30)

        @test lguys.calc_ρ.(prof, x ./ r_scale) * m_scale/r_scale^3 ≈ lguys.calc_ρ.(prof_scaled, x)
    end

end



@testset "ExpCusp" begin
    @testset "ρ" begin
        profile = lguys.ExpCusp(1, 1)

        @test lguys.calc_ρ(profile, 0) ≈ Inf
        @test lguys.calc_ρ(profile, 1) ≈ exp(-1) / 4π
        @test lguys.calc_ρ(profile, 2) ≈ exp(-2) / 2 / 4π
        @test lguys.calc_ρ(profile, Inf) ≈ 0
    end

    @testset "M" begin
        profile = lguys.ExpCusp(2.8, 1.11)
        x = 10 .^ LinRange(-1, 2, 10)

        M1 = lguys.calc_M.(profile, x)
        M2 = lguys.calc_M_from_ρ.(profile, x)

        @test M1 ≈ M2 rtol=1e-5
    end


    @testset "consistency" begin
        profile = lguys.ExpCusp(0.22, 0.31)

        test_r_h(profile)
        test_M_tot(profile)
        test_to_zero(profile)
    end
    
    @testset "scale" begin
        prof = lguys.ExpCusp(4.2, 12.3)
        r_scale = 1.3
        m_scale = π/2

        prof_scaled = lguys.scale(prof, r_scale, m_scale)

        x = 10 .^ LinRange(-1, 1, 30)

        @test lguys.calc_ρ.(prof, x ./ r_scale) * m_scale/r_scale^3 ≈ lguys.calc_ρ.(prof_scaled, x)
    end
end


@testset "Exp2D" begin
    @testset "Σ" begin
        profile = lguys.Exp2D(1, 1)

        @test lguys.calc_Σ(profile, 0) ≈ 1 / 2π 
        @test lguys.calc_Σ(profile, 1) ≈ exp(-1) / 2π
        @test lguys.calc_Σ(profile, 2) ≈ exp(-2) / 2π 
        @test lguys.calc_Σ(profile, Inf) ≈ 0
    end

    @testset "ρ" begin
        profile = lguys.Exp2D(0.25, √7)
        
        x = 10 .^ LinRange(0, 1, 10)

        ρ1 = lguys.calc_ρ.(profile, x)
        ρ2 = lguys.calc_ρ_from_Σ.(profile, x)

        @test ρ1 ≈ ρ2 rtol=1e-5
    end

    @testset "M" begin
        profile = lguys.Exp2D(0.989, 1.35)
        x = 10 .^ LinRange(-1, 2, 10)

        M1 = lguys.calc_M.(profile, x)
        M2 = lguys.calc_M_from_ρ.(profile, x)

        @test M1 ≈ M2 rtol=1e-5
    end


    @testset "consistency" begin
        M = 0.983
        r_s = 1.2
        profile = lguys.Exp2D(0.987, 3.3)

        test_R_h(profile)
        test_r_h(profile)
        test_M_tot(profile)
        test_to_zero(profile)
    end
    
    @testset "scale" begin
        prof = lguys.Exp2D(2.2, 0.7)
        r_scale = 1.3
        m_scale = π/2

        prof_scaled = lguys.scale(prof, r_scale, m_scale)

        x = 10 .^ LinRange(-1, 1, 30)

        @test lguys.calc_ρ.(prof, x ./ r_scale) * m_scale/r_scale^3 ≈ lguys.calc_ρ.(prof_scaled, x)
    end
end


@testset "Exp3D" begin
    @testset "ρ" begin
        profile = lguys.Exp3D(1, 1)

        @test lguys.calc_ρ(profile, 0) ≈ 1 / 8π 
        @test lguys.calc_ρ(profile, 1) ≈ exp(-1) / 8π
        @test lguys.calc_ρ(profile, 2) ≈ exp(-2) / 8π 
        @test lguys.calc_ρ(profile, Inf) ≈ 0
    end

    @testset "ρ" begin
        profile = lguys.Exp3D(0.25, √7)
        
        x = 10 .^ LinRange(0, 1, 10)

        ρ1 = lguys.calc_ρ.(profile, x)
        ρ2 = lguys.calc_ρ_from_Σ.(profile, x)

        @test ρ1 ≈ ρ2 rtol=1e-5
    end

    @testset "M" begin
        profile = lguys.Exp3D(0.989, 1.35)
        x = 10 .^ LinRange(-1, 2, 10)

        M1 = lguys.calc_M.(profile, x)
        M2 = lguys.calc_M_from_ρ.(profile, x)

        @test M1 ≈ M2 rtol=1e-5
    end


    @testset "consistency" begin
        profile = lguys.Exp3D(1.5, 0.95)

        test_R_h(profile)
        test_r_h(profile)
        test_M_tot(profile)
        test_to_zero(profile)
    end

    @testset "scale" begin
        prof = lguys.Exp3D(2.3, 0.8)
        r_scale = 1.4
        m_scale = π/3

        prof_scaled = lguys.scale(prof, r_scale, m_scale)

        x = 10 .^ LinRange(-1, 1, 30)

        @test lguys.calc_ρ.(prof, x ./ r_scale) * m_scale/r_scale^3 ≈ lguys.calc_ρ.(prof_scaled, x)
    end
end


@testset "log cusp 2D Σ" begin
    @testset "ρ" begin
        profile = lguys.LogCusp2D(1, 1)
        @test lguys.calc_ρ(profile, 0) ≈ Inf
        @test lguys.calc_ρ(profile, 1) ≈ exp(-1) / 4π
        @test lguys.calc_ρ(profile, 2) ≈ exp(-2) / 8π
    end

    @testset "consistency" begin
        profile = lguys.LogCusp2D(2.12, 3.53)
        test_R_h(profile)
        test_r_h(profile)
        test_M_tot(profile)
        test_to_zero(profile, r=10_000)
    end


    @testset "M" begin
        profile = lguys.LogCusp2D(0.987, 3.3)
        x = 10 .^ LinRange(-1, 2, 10)

        M1 = lguys.calc_M.(profile, x)
        M2 = lguys.calc_M_from_ρ.(profile, x)

        @test M1 ≈ M2 rtol=1e-5
    end


    @testset "Σ" begin
        profile = lguys.LogCusp2D(1, 1)

        x = 10 .^ LinRange(0, 1, 10)
        Σ1 = lguys.calc_Σ.(profile, x)
        Σ2 = lguys.calc_Σ_from_ρ.(profile, x)

        @test Σ1 ≈ Σ2 rtol=1e-5
    end

    @testset "scale" begin
        prof = lguys.LogCusp2D(3.5, 6.2)
        r_scale = 1.8
        m_scale = π^1.5

        prof_scaled = lguys.scale(prof, r_scale, m_scale)

        x = 10 .^ LinRange(-1, 1, 30)

        @test lguys.calc_ρ.(prof, x ./ r_scale) * m_scale/r_scale^3 ≈ lguys.calc_ρ.(prof_scaled, x)
    end
end


@testset "King" begin
    @testset "Σ" begin
        M = 1
        r_s = 1
        r_t = 2
        profile = lguys.KingProfile(k=1, R_s=r_s, R_t=r_t)

        @test lguys.calc_Σ(profile, 0) ≈ (1 - (1+(r_t/r_s)^2)^(-1/2))^2

        @test lguys.calc_Σ(profile, r_t) ≈ 0
        @test lguys.calc_Σ(profile, 1.2*r_t) ≈ 0
        @test lguys.calc_Σ(profile, Inf) ≈ 0
    end

    @testset "ρ" begin
        profile = lguys.KingProfile(k=0.666, R_s=2.3, R_t=9.5)
        x = 10 .^ LinRange(-1, 1, 10)

        ρ1 = lguys.calc_ρ.(profile, x)
        ρ2 = lguys.calc_ρ_from_Σ.(profile, x)
        @test ρ1 ≈ ρ2 rtol=1e-5
    end

    @testset "consistency" begin
        profile = lguys.KingProfile(M=1.2, R_s=2.3, R_t=9.5)
        test_R_h(profile, 1.2)
        test_r_h(profile, 1.2)
        test_M_tot(profile, 1.2)
        test_to_zero(profile)
    end


    @testset "scale" begin
        prof = lguys.KingProfile(M=2.2, R_s=1.9, R_t=3.3)
        r_scale = 1.7
        m_scale = π^2

        prof_scaled = lguys.scale(prof, r_scale, m_scale)

        x = 10 .^ LinRange(-1, 1, 30)

        @test lguys.calc_ρ.(prof, x ./ r_scale) * m_scale/r_scale^3 ≈ lguys.calc_ρ.(prof_scaled, x)
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
                    kwargs["k"] = rand()
                elseif profile == :NFW
                    kwargs = Dict()
                    kwargs["M_s"] = 10^randn()
                    kwargs["r_s"] = 10rand()
                    kwargs["c"] = 10rand()
                elseif profile == :TruncNFW
                    kwargs = Dict()
                    kwargs["M_s"] = 10^randn()
                    kwargs["r_s"] = 10rand()
                    kwargs["c"] = 10rand()

                    kwargs["r_t"] = 10rand()
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



@testset "calc_σv_star_mean" begin
    @test false broken=true
end

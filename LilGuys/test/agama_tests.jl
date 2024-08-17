

@testset "agama - my nfw" begin
    pot_agama = LilGuys.AgamaPotential(type="nfw")

    pot_me = LilGuys.NFW(M_s=1, r_s=1)

    r = 10 .^ LinRange(-1, 1, 10)

    x = LilGuys.calc_Φ.(pot_me, r)

    y = Float64[]

    println(pot_agama)
    GC.@preserve pot_agama for a in r
        println(pot_agama)
        println("r = ", a)
        phi = LilGuys.calc_Φ(pot_agama, [a, 0., 0.], 0.) 
        println("phi = ", phi)
        push!(y, phi)
    end


    @test x ≈ y
end

const I = [1 0 0; 0 1 0; 0 0 1]

function random_matrix()
    return lguys.Rz_mat(rand() * 2π) * lguys.Ry_mat(rand() * 2π) * lguys.Rx_mat(rand() * 2π)
end

@testset "angular_distance" begin
    @testset "angular_distance simple cases" begin
        @test lguys.angular_distance(0, 0, 0, 0) ≈ 0
        @test lguys.angular_distance(0, 0, 0, 90) ≈ 90
        @test lguys.angular_distance(0, 0, 90, 0) ≈ 90
        @test lguys.angular_distance(0, 0, -90, 0) ≈ 90
        @test lguys.angular_distance(0, 0, 0, -90) ≈ 90
        @test lguys.angular_distance(0, 0, 180, 0) ≈ 180
        @test lguys.angular_distance(0, 0, 0, 180) ≈ 180
    end


    @testset "angular_distance identity" begin
        # tricy almost same test
        @test lguys.angular_distance.([14.479505360741172], 
                                      [-33.99473291943769], [14.479505360741172], [-33.99473291943769]) == [0.]

        N = 1000
        ra = rand(N) * 360
        dec = rand(N) * 180 .- 90

        @test lguys.angular_distance.(ra, dec, ra, dec) ≈ zeros(N)  atol=1e-4
    end

    @testset "random" begin
        # use rotations to setup random two points with know distance
        for _ in 1:10
            dist = rand() * 180

            mat = random_matrix()

            alpha1, beta1 = lguys.rotate_sky(rand() * 360, 90-dist, mat)
            alpha2, beta2 = lguys.rotate_sky(0, 90, mat)

            @test lguys.angular_distance(alpha1, beta1, alpha2, beta2) ≈ dist
            @test lguys.angular_distance(alpha2, beta2, alpha1, beta1) ≈ dist
            @test lguys.angular_distance(alpha1, beta1, alpha1, beta1) ≈ 0
        end
    end

    @testset "random, small angle" begin
        for _ in 1:10
            dist = rand() / 3600 # 1 arcmin

            mat = random_matrix()

            @test lguys.angular_distance(0, 90, 0, 90-dist) ≈ dist atol=1e-10

            alpha1, beta1 = lguys.rotate_sky(rand() * 360, 90-dist, mat)
            alpha2, beta2 = lguys.rotate_sky(0, 90, mat)

            @test lguys.angular_distance(alpha1, beta1, alpha2, beta2) ≈ dist atol=1e-10
            @test lguys.angular_distance(alpha2, beta2, alpha1, beta1) ≈ dist atol=1e-10
        end
    end


    @testset "random, very small angle" begin
        for _ in 1:10
            dist = 1e-8 * rand() # 1 arcmin

            mat = random_matrix()

            @test lguys.angular_distance(0, 90, 0, 90-dist) ≈ dist atol=1e-10

            alpha1, beta1 = lguys.rotate_sky(rand() * 360, 90-dist, mat)
            alpha2, beta2 = lguys.rotate_sky(0, 90, mat)

            @test lguys.angular_distance(alpha1, beta1, alpha2, beta2) ≈ dist atol=1e-10
            @test lguys.angular_distance(alpha2, beta2, alpha1, beta1) ≈ dist atol=1e-10
        end
    end
end



@testset "to_tangent" begin
    @testset "special cases" begin
        α_0 = 0
        δ_0 = 0

        α = Float64[1, 0, -1, 0]
        δ = Float64[0, 1, 0, -1]

        xi, eta = lguys.to_tangent(α, δ, α_0, δ_0)

        @test xi ≈ α atol=1e-2
        @test eta ≈ δ atol=1e-2
    end

    @testset "random" begin
        for _ in 1:100
            α_0 = 0
            δ_0 = 0

            # want to preseve north
            mat = lguys.Rz_mat(rand() * 2π) * lguys.Ry_mat(π - 2π * rand()) 
            α_0, δ_0 = lguys.rotate_sky(α_0, δ_0, mat)



            # identity
            α = 0
            δ = 0
            α_p, δ_p = lguys.rotate_sky(α, δ, mat)
            xi, eta =  lguys.to_tangent(α_p, δ_p, α_0, δ_0)
            @test xi ≈ 0
            @test eta ≈ 0

            #NaN
            α = [90 + rand() * 180, 0]
            δ = [0, 90 + rand() * 180]

            α_p, δ_p = lguys.rotate_sky(α, δ, mat)
            xi, eta =  lguys.to_tangent(α_p, δ_p, α_0, δ_0)
            @test all(isnan.(xi))
            @test all(isnan.(eta))


            # eta hat
            α = randn()
            δ = 0
            α_p, δ_p = lguys.rotate_sky(α, δ, mat)
            xi, eta =  lguys.to_tangent(α_p, δ_p, α_0, δ_0)

            #transformations complicate the direction
            xi_sign = mod.(α_p - α_0, 360) .> 180 ? -1 : 1
            @test xi ≈ abs(α) * xi_sign atol=1e-2
            @test eta ≈ 0 atol=1e-12

            # xi hat
            α = 0
            δ = randn()
            α_p, δ_p = lguys.rotate_sky(α, δ, mat)
            xi, eta =  lguys.to_tangent(α_p, δ_p, α_0, δ_0)
            @test xi ≈ 0 atol=1e-12
            @test eta ≈ abs(δ) * sign(δ_p - δ_0) atol=1e-2

            # magnitude
            # TODO:
        end
    end

end


@testset "rotate_sky" begin
    @testset "simple cases" begin
        ra = [0, 90, 180, 270]
        dec = [0, 0, 0, 0]

        ra_p, dec_p = lguys.rotate_sky(ra, dec, lguys.Rz_mat(π/2))

        @test mod.(ra_p, 360) ≈ [90, 180, 270, 0]
        @test mod.(dec_p, 360) ≈ [0, 0, 0, 0]
    end

    @testset "identity" begin
        for _ in 1:100
            ra = rand() * 360
            dec = rand() * 180 .- 90

            @test lguys.angular_distance(lguys.rotate_sky(ra, dec, I)..., ra, dec) ≈ 0 atol=1e-12
        end
    
    end
end


@testset "unit_vector" begin
    @testset "basis vectors" begin
        @test lguys.unit_vector(0, 0) ≈ [1, 0, 0]
        @test lguys.unit_vector(90, 0) ≈ [0, 1, 0]
        @test lguys.unit_vector(180, 0) ≈ [-1, 0, 0]
        @test lguys.unit_vector(270, 0) ≈ [0, -1, 0]
        @test lguys.unit_vector(0, 90) ≈ [0, 0, 1]
        @test lguys.unit_vector(0, -90) ≈ [0, 0, -1]
    end

    @testset "simple cases" begin
        @test lguys.unit_vector(45, -60) ≈ [√2/4, √2/4, -√3/2]
    end

    @testset "normalization" begin
        for _ in 1:100
            ra = rand() * 360
            dec = rand() * 180 .- 90

            v = lguys.unit_vector(ra, dec)

            # normalization
            @test lguys.norm(v) ≈ 1

            # periodicity
            v2 = lguys.unit_vector(ra + 360, dec)
            @test v ≈ v2
            v2 = lguys.unit_vector(ra, dec + 360)
            @test v ≈ v2
        end
    end

end


@testset "cartesian_to_sky" begin
    @test false broken = true
end


@testset "Rx_mat" begin
    # for 90 degrees
    # yhat -> zhat 
    # zhat -> -yhat
    H = lguys.Rx_mat(π/2)
    @test H * [0, 1, 0] ≈ [0, 0, 1]
    @test H * [0, 0, 1] ≈ [0, -1, 0]
    @test H * [1, 0, 0] ≈ [1, 0, 0]

    H = lguys.Rx_mat(π)
    @test H * [0, 1, 0] ≈ [0, -1, 0]
    @test H * [0, 0, 1] ≈ [0, 0, -1]
    @test H * [1, 0, 0] ≈ [1, 0, 0]

    H = lguys.Rx_mat(π/6)
    @test H * [0, 1, 0] ≈ [0, √3/2, 1/2]
    @test H * [0, 0, 1] ≈ [0, -1/2, √3/2]
    @test H * [1, 0, 0] ≈ [1, 0, 0]
end


@testset "Ry_mat" begin
    # for 90 degrees
    # xhat -> zhat 
    # zhat -> -xhat
    H = lguys.Ry_mat(π/2)
    @test H * [1, 0, 0] ≈ [0, 0, -1]
    @test H * [0, 0, 1] ≈ [1, 0, 0]
    @test H * [0, 1, 0] ≈ [0, 1, 0]

    H = lguys.Ry_mat(π)
    @test H * [1, 0, 0] ≈ [-1, 0, 0]
    @test H * [0, 0, 1] ≈ [0, 0, -1]
    @test H * [0, 1, 0] ≈ [0, 1, 0]

    H = lguys.Ry_mat(π/6)
    @test H * [1, 0, 0] ≈ [√3/2, 0, -1/2]
    @test H * [0, 0, 1] ≈ [1/2, 0, √3/2]
    @test H * [0, 1, 0] ≈ [0, 1, 0]
end



@testset "Rz_mat" begin
    # for 90 degrees
    # xhat -> yhat 
    # yhat -> -xhat
    H = lguys.Rz_mat(π/2)
    @test H * [1, 0, 0] ≈ [0, 1, 0]
    @test H * [0, 1, 0] ≈ [-1, 0, 0]
    @test H * [0, 0, 1] ≈ [0, 0, 1]

    H = lguys.Rz_mat(π)
    @test H * [1, 0, 0] ≈ [-1, 0, 0]
    @test H * [0, 1, 0] ≈ [0, -1, 0]
    @test H * [0, 0, 1] ≈ [0, 0, 1]

    H = lguys.Rz_mat(π/6)
    @test H * [1, 0, 0] ≈ [√3/2, 1/2, 0]
    @test H * [0, 1, 0] ≈ [-1/2, √3/2, 0]
    @test H * [0, 0, 1] ≈ [0, 0, 1]
end



@testset "rotation identity" begin

    for θ in [-2π, 0, 2π]
        @test lguys.Rx_mat(θ) ≈ I
        @test lguys.Ry_mat(θ) ≈ I
        @test lguys.Rz_mat(θ) ≈ I
    end
end


@testset "rotation inverse" begin
    for _ in 1:20
        θ = rand() * 2π

        for mat in [lguys.Rx_mat, lguys.Ry_mat, lguys.Rz_mat]
            H = mat(θ)
            H_inv = mat(-θ)

            @test H * H_inv ≈ I
        end
    end
end


@testset "rotation periodic" begin
    for _ in 1:20
        θ = rand() * 2π

        for mat in [lguys.Rx_mat, lguys.Ry_mat, lguys.Rz_mat]
            H = mat(θ)
            H_inv = mat(θ + 2π)

            @test H ≈ H_inv
        end
    end
end

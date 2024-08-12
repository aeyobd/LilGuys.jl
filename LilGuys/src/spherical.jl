import LinearAlgebra: det

"""
    to_tangent(α, δ, α_0, δ_0)

Computes the tangent plane coordinates of a point (α, δ) with respect to a reference point (α_0, δ_0).
α and δ may be reals or vectors.
Both input and output angles are assumed to be in degrees.
Returns NaN for points with angular distances greater than 90 degrees.
"""
function to_tangent(α::Real, δ::Real, α_0::Real, δ_0::Real)
    xi, eta = _to_tangent(α, δ, α_0, δ_0)
    if angular_distance(α, δ, α_0, δ_0) > 90
        xi = NaN
        eta = NaN
    end
    return xi, eta
end



function to_tangent(α::AbstractVector, δ::AbstractVector, α_0::Real, δ_0::Real)
    xi = similar(α)
    eta = similar(δ)

    for i in eachindex(α)
        xi[i], eta[i] = to_tangent(α[i], δ[i], α_0, δ_0)
    end

    return xi, eta
end



function _to_tangent(α, δ, α_0, δ_0)
    denom = sind(δ) * sind(δ_0) + cosd(δ) * cosd(δ_0) * cosd(α - α_0)

    xi_num = cosd(δ) * sind(α - α_0)
    eta_num = sind(δ) * cosd(δ_0) - cosd(δ) * sind(δ_0) * cosd(α-α_0) 

    xi = rad2deg(xi_num/denom)
    eta = rad2deg(eta_num / denom)
    return xi, eta
end


"""
    angular_distance(α1, δ1, α2, δ2)

Computes the angular distance in degrees between two points on the sky, 
assuming RA/DEC and in degrees.
Uses Haversine formula
"""
function angular_distance(α1::Real, δ1::Real, α2::Real, δ2::Real)
    dα  = α2 - α1
    dδ = δ2 - δ1

    a = sind(dδ/2)^2 + cosd(δ1) * cosd(δ2) * sind(dα/2)^2
    c = 2 * atan(sqrt(a), sqrt(1-a))
    return c * 180 / π
end



"""
    unit_vector(ra, dec)

Returns the unit vector(s) pointing at the position(s) (ra, dec) in degrees on the sky. ra and dec may be reals or vectors, but should be in degrees.
"""
function unit_vector(ra::Real, dec::Real)
    x = cosd(dec) * cosd(ra)
    y = cosd(dec) * sind(ra)
    z = sind(dec)
    return [x, y, z]
end


function unit_vector(ra::Vector{<:Real}, dec::Vector{<:Real}) 
    x = cosd.(dec) .* cosd.(ra)
    y = cosd.(dec) .* sind.(ra)
    z = sind.(dec)
    return [x y z]'
end


"""
    cartesian_to_sky(x, y, z)

Converts cartesian coordinates to spherical coordinates (RA, DEC, r)
"""
function cartesian_to_sky(x::T, y::T, z::T) where T <: Union{Real, AbstractVector}
    R = sqrt.(x.^2 + y.^2)
    r = sqrt.(x.^2 + y.^2 + z.^2)

    ra = mod.(atand.(y, x), 360) # atan2 -> 0 to 360
    dec = atand.(z ./ R) # atan -> -90 to 90

    return ra, dec, r
end


function cartesian_to_sky(mat::AbstractMatrix)
    @assert size(mat, 1) == 3
    return cartesian_to_sky(mat[1, :], mat[2, :], mat[3, :])
end

function cartesian_to_sky(v::AbstractVector)
    @assert length(v) == 3
    return cartesian_to_sky(v[1], v[2], v[3])
end

"""
    sky_to_cartesian(ra, dec, r=1)

Converts spherical coordinates to cartesian coordinates
"""
function sky_to_cartesian(ra::Real, dec::Real, r=1)
    return r .* unit_vector(ra, dec)
end


"""
    rotate_sky(ra, dec, R)

Given a point (ra, dec) on the sky, rotates it by the rotation matrix R
"""
function rotate_sky(ra, dec, R)
    @assert size(R) == (3, 3)
    @assert det(R) ≈ 1

    x_mat = unit_vector(ra, dec)
    x_mat = R * x_mat

    ra_p, dec_p, _ = cartesian_to_sky(x_mat)

    return ra_p, dec_p
end



"""
    Rx_mat(u)

rotation matrix around x axis (radians)
"""
function Rx_mat(u::Real)
    c = cos(u)
    s = sin(u)
    return [1  0  0
            0  c -s
            0  s  c]
end



"""
    Ry_mat(v)

rotation matrix around y axis (radians)
"""
function Ry_mat(v::Real)
    c = cos(v)
    s = sin(v)
    return [c  0  s
            0  1  0
           -s  0  c]
end


"""
    Rz_mat(w)

rotation matrix around z axis (radians)
"""
function Rz_mat(w::Real)
    c = cos(w)
    s = sin(w)

    return [c -s  0
            s  c  0
            0  0  1]
end



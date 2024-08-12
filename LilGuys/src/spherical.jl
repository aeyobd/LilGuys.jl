"""
    to_tangent(α, δ, α_0, δ_0)

Computes the tangent plane coordinates of a point (α, δ) with respect to a reference point (α_0, δ_0).
α and δ may be reals or vectors
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
    eta_num = sind(δ_0) * cosd(δ) * cosd(α-α_0) - cosd(δ_0) * sind(δ)
    xi_num = cosd(δ) * sind(α - α_0)
    xi = rad2deg(xi_num/denom)
    eta = -rad2deg(eta_num / denom)
    return xi, eta
end


"""
    angular_distance(α1, δ1, α2, δ2)

Computes the angular distance in degrees between two points on the sky, 
assuming RA/DEC and in degrees.
"""
function angular_distance(α1::Real, δ1::Real, α2::Real, δ2::Real)
    a = sind(δ1) * sind(δ2) + cosd(δ1) * cosd(δ2) * cosd(α1 - α2)
    if abs(a) > nextfloat(1.0, 5)
        raise(DomainError("a greater than 1: $a"))
    end

    if abs(a) > 1
        a = sign(a) * 1
    end

    return acosd(a)
end



"""
    unit_vector(ra, dec)

Returns the unit vector(s) pointing at the position(s) (ra, dec) in degrees on the sky. ra and dec may be reals or vectors, but should be in degrees.
"""
function unit_vector(ra::F, dec::F)
    x = cosd(dec) * cosd(ra)
    y = cosd(dec) * sind(ra)
    z = sind(dec)
    return [x, y, z]
end


function unit_vector(ra::Vector{F}, dec::Vector{F}) where F
    x = cosd.(dec) .* cosd.(ra)
    y = cosd.(dec) .* sind.(ra)
    z = sind.(dec)
    return [x y z]'
end


"""
    cartesian_to_sky(x, y, z)

Converts cartesian coordinates to spherical coordinates (RA, DEC, r)
"""
function cartesian_to_sky(x, y, z)

    R = sqrt(x^2 + y^2)
    r = sqrt(x^2 + y^2 + z^2)

    ra = mod(atand(y, x), 360) # atan2 -> 0 to 360
    dec = atand(z, R) # atan -> -90 to 90

    return [ra, dec, r]
end



function cartesian_to_sky(x::AbstractVector, y::AbstractVector, z::AbstractVector)
    R = sqrt.(x.^2 + y.^2)
    r = sqrt.(x.^2 + y.^2 + z.^2)

    ra = mod.(atand.(y, x), 360) # atan2 -> 0 to 360
    dec = atand.(z ./ R) # atan -> -90 to 90

    return [ra dec r]
end


function sky_to_cartesian(ra, dec, r)
    return r .* unit_vector(ra, dec)
end


"""
    Rx_mat(u)

rotation matrix around x axis (radians)
"""
function Rx_mat(u::Real)
    c = cos(u)
    s = sin(u)
    return [1  0  0
            0  c  s
            0 -s  c]
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


"""
Rotation matrix 

R = Rz(w) * Ry(v) * Rx(u)
"""
function R_mat(u::Real, v::Real, w::Real)
    return Rz_mat(w) * Ry_mat(v) * Rx_mat(u)
end



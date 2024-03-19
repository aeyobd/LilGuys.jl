import NearestNeighbors as nn


struct KDE
    x::Matrix{Float64}
    h::Vector{Float64}
    kernel::Function
    dim::Int
    trunc::Float64 = 5
    _tree
end


function (kde::KDE)(x::Vector{Float64})

end

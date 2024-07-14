using Glob



Base.@kwdef struct Output <: AbstractArray{Snapshot, 1}
    """ The hdf5 file containing the snapshots """
    h5file::HDF5.File

    """ The times of the snapshots (code units) """
    times::Vector{F}

    """The index of the snapshots"""
    index::Vector{String}

    """ The (calculated) position centres of the snapshots """
    x_cen::Matrix{F} = zeros(F, 3, length(times))

    """ The (calculated) velocity centres of the snapshots """
    v_cen::Matrix{F} = zeros(F, 3, length(times))

    """ The stellar weights of the system """
    weights::OptVector = nothing

    """ The gravitational softening length """
    softening::F = NaN
end


function Base.finalize(out::Output)
    close(out.h5file)
end


""" Find the output file in the given path."""
function _find_output_filename(path::String)
    if splitext(basename(path))[2] == ".hdf5"
        return path
    elseif isdir(path)
        if "combined.hdf5" ∈ readdir(path)
            return joinpath(path, "combined.hdf5")
        elseif ("out" ∈ readdir(path)) && ("combined.hdf5" ∈ readdir(joinpath(path, "out")))
            return joinpath(path, "out", "combined.hdf5")
        else
            error("combined.hdf5 not found in $(path)")
        end
    end
end


"""
Create an output object from the given filename.
The stellar weights can be provided as an optional argument.
"""
function Output(filename::String; weights=nothing)

    filename = _find_output_filename(filename)

    file = h5open(filename, "r")
    names = keys(file)
    snap_names = filter(x -> startswith(x, "snap"), names)

    Nt = length(snap_names)
    index = Vector{String}(undef, Nt)
    times = Vector{F}(undef, Nt)

    for i in 1:Nt
        index[i] = "snap$(i-1)"
        if index[i] ∉ names
            error("$(index[i]) not found in file")
        end
        header = get_header(file[index[i]])
        times[i] = header["Time"]
    end

    if "x_cen" in names
        x_cen = file["x_cen"][:, :]
    else
        x_cen = zeros(F, 3, Nt)
    end

    if "v_cen" in names
        v_cen = file["v_cen"][:, :]
    else
        v_cen = zeros(F, 3, Nt)
    end

    out = Output(;h5file=file, times=times, index=index, x_cen=x_cen, v_cen=v_cen, weights=weights)
    return out
end


function Base.size(out::Output)
    N = length(out.index)
    return (N,)
end

Base.IndexStyle(::Type{<:Output}) = IndexCartesian()


function Base.getindex(out::Output, i::Int)
    snap =  Snapshot(out.h5file[out.index[i]])
    snap.x_cen = out.x_cen[:, i]
    snap.v_cen = out.v_cen[:, i]

    if !isnothing(out.weights)
        if length(snap.index) !== length(out.weights)
            println("warning: Weights and snapshot index have different lengths: $(length(snap.index)) vs $(length(out.weights))")

        else
            if isperm(snap.index)
                idx = snap.index
            else
                idx = invperm(sortperm(snap.index))
            end

            snap.weights = out.weights[idx]
        end
    end

    return snap
end



"""
    extract(snap::Snapshot, symbol::Symbol, idx::Int)

Extract the value of a field from a snapshot (at a given index). Returns a list as sorted by index
"""
function extract(snap::Snapshot, symbol, idx=(:))
    idx_sort = sortperm(snap.index)
    attr = getfield(snap, symbol)
    val = attr[idx_sort[idx]]
    return val
end


"""
    extract_vector(snap::Snapshot, symbol::Symbol, idx::Int)

Extract the value of a field from a snapshot (at a given index). Returns a list as sorted by index
"""
function extract_vector(snap::Snapshot, symbol, idx=(:))
    idx_sort = sortperm(snap.index)
    attr = getfield(snap, symbol)
    val = attr[:, idx_sort[idx]]
    return val
end




"""
Extracts the given symbol from the output at the given index
"""
function extract(out::Output, symbol::Symbol, idx::Int)
    Nt = length(out)
    result = Array{F}(undef, Nt)

    for i in 1:Nt
        h5f = out.h5file[out.index[i]]
        snap_idx = get_vector(h5f, h5vectors[:index])
        idx_sort = sortperm(snap_idx)
        result[i] = h5f["PartType1/$(h5vectors[symbol])"][idx_sort[idx]]
    end
    return result
end




"""
Extracts the given symbol from the output at the given index
"""
function extract_vector(out::Output, symbol::Symbol, idx::Int; dim::Int=3)
    Nt = length(out)
    result = Array{F}(undef, dim, Nt)

    for i in 1:Nt
        h5f = out.h5file[out.index[i]]
        snap_idx = get_vector(h5f, h5vectors[:index])
        idx_sort = sortperm(snap_idx)
        result[:, i] = h5f["PartType1/$(h5vectors[symbol])"][:, idx_sort[idx]]
    end
    return result
end



function extract(out::Output, symbol::Symbol, idx=(:))
    if idx == (:)
        idx = 1:length(out[1].index)
    elseif idx isa BitArray
        idx = findall(idx)
    end
    Np = length(idx)
    Nt = length(out)
    result = Array{F}(undef, Np, Nt)

    for i in 1:Nt
        snap = out[i]
        result[:, i] .= extract(snap, symbol, idx)
    end

    return result
end


function extract_vector(out::Output, symbol::Symbol, idx=(:))
    if idx == (:)
        idx = 1:length(out[1].index)
    elseif idx isa BitArray
        idx = findall(idx)
    end
    Np = length(idx)
    Nt = length(out)
    result = Array{F}(undef, 3, Np, Nt)

    for i in 1:Nt
        snap = out[i]
        result[:, :, i] .= extract_vector(snap, symbol, idx)
    end

    return result
end




function peris_apos(out::Output; verbose::Bool=false)
    r0 = calc_r(out[1].positions)
    peris = apos = r0

    idx0 = sort(out[1].index)

    if verbose
        println("begining peri apo calculation")
    end

    N = length(out)
    for i in 2:N
        if verbose && i % 100 == 0
            print("$(i)/$(N)\r")
        end
        
        snap = out[i]
        idx = sortperm(snap.index)

        r = calc_r(snap.positions[:, idx])
        apos = max.(apos, r)
        peris = min.(peris, r)

    end
    
    if verbose
        println("completed peri apo calculation")
    end
    return idx0, peris, apos  
end



function Base.show(io::IO, out::Output)
    print(io, "<output with $(length(out)) snapshots of $(length(out[1])) particles>")
    return io
end


function Base.show(io::IO, mime::MIME"text/plain", out::Output)
    print(io, out)
    return io
end

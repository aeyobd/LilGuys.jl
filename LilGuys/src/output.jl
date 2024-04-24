using Glob



Base.@kwdef struct Output <: AbstractArray{Snapshot, 1}
    h5file::HDF5.File
    times::Vector{F}
    softening::F
    index::Vector{String}
end


function Base.finalize(out::Output)
    close(out.h5file)
end



function Output(filename::String)
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

    
    out = Output(file, times, NaN, index)
    return out
end


function Base.size(out::Output)
    N = length(out.index)
    return (N,)
end

Base.IndexStyle(::Type{<:Output}) = IndexCartesian()

function Base.getindex(out::Output, i::Int)
    return Snapshot(out.h5file[out.index[i]])
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

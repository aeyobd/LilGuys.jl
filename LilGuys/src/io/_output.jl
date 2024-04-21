using Glob


Base.@kwdef struct Output <: AbstractArray{Snapshot, 1}
    snapshots::Vector{Snapshot}
    times::Vector{F}
    softening::F
end


Base.@kwdef struct Output2 <: AbstractArray{Snapshot, 1}
    h5file::HDF5File
    times::Vector{F}
    softening::F
end




function Output(directory::String)

    filenames = glob("snapshot*.hdf5", directory)

    snap0 = Snapshot(filenames[1], mmap=true)
    idx0 = sort(snap0.index)
    Np = length(idx0)
    Nt = length(filenames)

    snapshots = Snapshot[]
    time = F[]
    softening = get_epsilon(directory)

    for i in 1:Nt
        filename = filenames[i]
        snap = Snapshot(filename, mmap=true)
        push!(snapshots, snap)

        t = snap.header["Time"]
        push!(time, t)

        idx = snap.index
        if sort(idx) != idx0
            error("Non-constant index")
        end
    end
    time_index = sortperm(time)

    
    out = Output(snapshots[time_index], time[time_index], softening)
    return out
end

function out_from_file(filename::String)
    snap = Snapshot(filename, mmap=true)
    out = Output([snap], [snap.header["Time"]], get_epsilon(dirname(filename)))
    return out
end


function Base.size(out::Output)
    N = length(out.snapshots)
    return (N,)
end


Base.IndexStyle(::Type{<:Output}) = IndexCartesian()

function Base.getindex(out::Output, i::Int)
    return out.snapshots[i]
end



function get_epsilon(dir::String)
    filename = joinpath(dir, "parameters-usedvalues")
    for line in eachline(filename)
        if startswith(line, "SofteningComovingClass0")
            m = collect(eachmatch(r"\d*\.?\d+", line))[end]
            return parse(F, m.match)
        end
    end
    return NaN
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
        snap = out[i]
        result[i] = extract(snap, symbol, idx)
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
        snap = out[i]
        result[:, i] = extract_vector(snap, symbol, idx)
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

        if any(snap.index[idx] .!= idx0)
            error("Non-constant index")
        end
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

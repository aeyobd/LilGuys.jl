using Glob


Base.@kwdef struct Output <: AbstractArray{Snapshot, 1}
    snapshots::Vector{Snapshot}
    times::Vector{F}
    softening::F
    particle_index::Matrix{Int}
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
    particle_index = Matrix{Int}(undef, Nt, Np)

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

        particle_index[i, :] = sortperm(idx)
    end
    time_index = sortperm(time)

    
    out = Output(snapshots[time_index], time[time_index], softening, particle_index[time_index, :])
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

function extract(out::Output, symbol::Symbol, idx::Int)
    Nt = length(out)
    if symbol ∈ (:positions, :velocities, :accelerations)
        var_is_vec = true
        result = Array{F}(undef, 3, Nt)
    else
        result = Array{F}(undef, Nt)
        var_is_vec = false
    end

    for i in 1:Nt
        snap = out[i]
        idx_sort = sortperm(snap.index)[idx]
        attr = getfield(snap, symbol)
        if var_is_vec
            val = attr[:, idx_sort]
            result[:, i] .= val
        else
            val = attr[idx_sort]
            result[i] = val
        end
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
    if symbol ∈ (:positions, :velocities, :accelerations)
        var_is_vec = true
        result = Array{F}(undef, 3, Np, Nt)
    else
        result = Array{F}(undef, Np, Nt)
        var_is_vec = false
    end

    for i in 1:Nt
        snap = out[i]
        idx_sort = sortperm(snap.index)[idx]
        attr = getfield(snap, symbol)
        if var_is_vec
            val = attr[:, idx_sort]
            result[:, :, i] .= val
        else
            val = attr[idx_sort]
            result[:, i] .= val
        end
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

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
    Nt = length(out.snapshots)
    Np = length(out.snapshots[1])
    return (Nt, Np)
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
    return nothing
end


function extract(out::Output, func::Function, idx_P=(:))
    Nt, Np = size(out)
    return [func(out[i])]
end


function peris_apos(out::Output)
    Nt, Np = size(out)
    
    r0 = [norm(p.pos) for p in out[1, :]]
    peris = r0
    apos = r0

    println("begining")
    for i in 2:Nt
        if i % 100 == 0
            print("$(i)/$(Nt)\r")
        end
        
        ps = out[i, :]
        r = [norm(p.pos) for p in ps]
        apos = max.(apos, r)
        peris = min.(peris, r)
    end
    
    println("\r finished")
    return peris, apos  
end

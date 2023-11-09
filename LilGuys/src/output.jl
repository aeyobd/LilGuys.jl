using Glob


Base.@kwdef struct Output <: AbstractArray{Particle, 2}
    snapshots::Vector{Snapshot}
    time::Vector{F}
    softening::F
end

function Output(directory::String)
    snapshots = Snapshot[]
    time = F[]
    path = joinpath(directory, "*.hdf5")
    softening = get_epsilon(directory)

    out = Output(snapshots, time, softening)
    for filename in glob(path)
        snap = Snapshot(filename, mmap=true)
        push!(out.snapshots, snap)
        t = snap.header["Time"]
        push!(out.time, t)
    end

    return out
end

function Base.size(out::Output)
    Nt = length(out.snapshots)
    Np = length(out.snapshots[1])
    return (Nt, Np)
end


Base.IndexStyle(::Type{<:Output}) = IndexCartesian()

function Base.getindex(out::Output, i, j)
    return out.snapshots[i][j]
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




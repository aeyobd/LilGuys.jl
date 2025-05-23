import DataFrames: DataFrame
import HDF5
using OrderedCollections

h5open = HDF5.h5open


"""
    read_hdf5_table(filename, path="/")

Load an HDF5 file and return a DataFrame using the specified path.
"""
function read_hdf5_table(filename::String, path="/")
    df = DataFrame()

    h5open(filename, "r") do f
        for k in keys(f[path])
            if k == "Header"
                continue
            end

            val = HDF5.read(f[path * "/" * k])

            df[!, k] = val
        end
    end

    return df
end



"""
    write_hdf5_table(filename, dataframe; path="/", overwrite=false)

Write a DataFrame to an HDF5 file.
"""
function write_hdf5_table(filename::String, frame::DataFrame; path="", overwrite=false)
    if isfile(filename) && !overwrite
        throw(ArgumentError("File already exists: $filename"))
    elseif isfile(filename) && overwrite
        rm(filename, force=true)
    end

    for col in names(frame)
        val = frame[:, col]

        if val isa BitVector
            val = Int.(val)
        end

        HDF5.h5write(filename, path * "/" * col, val)
    end
end


"""
Converts a Dataframe to a Dict{String, Any} object.
"""
function to_dict(frame::DataFrame)
    df = OrderedDict(String(name) => frame[:, name] for name in names(frame))

    return df
end


"""
    remove_nonscalars(d)

return a copy of the dictionary except only retaining
strings, floats, and splatting arrays
"""
function remove_nonscalars(d::AbstractDict)
    d_clean = copy(d)
    for (k, val) in d
        if val isa AbstractDict
            @debug "removing dict value $k"
            delete!(d_clean, k)
        elseif val isa AbstractArray
            if eltype(val) <: Union{Real, String}
                for (i, v) in enumerate(val)
                    d_clean[Symbol("$(k)_$i")] = v
                end
                delete!(d_clean, k)
            else
                @warn "type $(typeof(val)) not implemented"
            end
        elseif val isa String
        elseif val isa Real
        else
            @warn "type $(typeof(val)) not implemented"
        end
    end
    
    return d_clean
end


"""
    to_frame(structs)

Takes a vector of scalar structs and converts them to a dataframe
"""
function to_frame(A::AbstractVector{<:T}) where T
    cols = propertynames(A[1])
    df = DataFrame()
    for col in cols
        val = getproperty.(A, col)
        
        if val[1] isa Union{Real, String}
            df[!, Symbol(col)] = val
        elseif val[1] isa AbstractVector{<:Real}
            @assert all(x->length(x) == length(val[1]), val)
            for i in eachindex(val[1])
                colname = Symbol("$(col)_$i")
                @assert colname ∉ names(df)
                df[!, colname] = [v[i] for v in val]
            end
        else
            @info "Skipping column $col"
        end 
    end

    return df
end


"""
    to_structs(df, type)

Converts a dictionary into a list of structs of the type
assuming the struct can be initialized with kwargs equal
to the names of the dataframe.
"""
function to_structs(df::DataFrame, T::Type)
    out = Vector{T}(undef, size(df, 1))

    for (i, row) in enumerate(eachrow(df))
        kwargs = Dict(Symbol(k) => row[k] for k in names(df))
        out[i] = T(; kwargs...)
    end

    return out
end


"""
    set_header!(h5f, header)

Sets the header of an HDF5 file using the given dictionary.
"""
function set_header!(h5f::HDF5.H5DataStore, header::Dict{String,<:Any})
    if "Header" ∉ keys(h5f)
        HDF5.create_group(h5f, "Header")
    end
    h5_header = h5f["Header"]
    for (key, val) in header
        set_attr!(h5_header, key, val)
    end
end


"""
    set_attrs!(h5f, header)

Sets the attributes of the given hdf5 data store to the provided dict
"""
function set_attrs!(h5f::HDF5.H5DataStore, header::Dict{String, <:Any})
    for (key, val) in header
        set_attr!(h5f, key, val)
    end
end
    


"""
    set_header_attr(h5f, key, val)

Sets an attribute in the header of an HDF5 file.
"""
function set_attr!(h5f::HDF5.H5DataStore, key::String, val)
    header = HDF5.attrs(h5f)
    header[key] = val
end



"""
    get_header(h5f)

Returns the header of an HDF5 file as a dictionary.
"""
function get_header(h5f::HDF5.H5DataStore)
    return get_attrs(h5f["Header"])
end


"""
    get_attrs(h5f, group)
"""
function get_attrs(h5f::HDF5.H5DataStore)
    return OrderedDict(HDF5.attrs(h5f))
end


"""
    get_vector(h5f, key; mmap=false, group="PartType1")

Returns a vector from an HDF5 file.
"""
function get_vector(h5f::HDF5.H5DataStore, key::String; mmap=false)
    if mmap
        return HDF5.readmmap(h5f[key])
    else
        return read(h5f[key])
    end
end



"""
    set_vector!(h5f, key, val; group="PartType1")

Sets a vector in an HDF5 file.
"""
function set_vector!(h5f::HDF5.File, key::String, val)
    h5f[key] = val
end



"""
    write_structs_to_hdf5
"""
function write_structs_to_hdf5(filename::String, structs::AbstractVector{<:Pair{String, <:Any}}; 
        overwrite=true)
    if overwrite && isfile(filename)
        rm(filename, force=true)
    end

    h5open(filename, "w") do f
        for (label, obj) in structs
            write_struct_to_hdf5(f, obj, group=label)
        end
    end

end



"""
    read_structs_from_hdf5(filename, T)

Reads a vector of structs from an HDF5 file. 
Expects each struct to be stored in a separate group and
T to be initializable by kwdef by the combination of 
values stored in each group.
""" 
function read_structs_from_hdf5(filename::String, T; use_measurements=:auto)
    structs = Pair{String, T}[]

    h5open(filename, "r") do f
        for k in keys(f)
            @debug "reading key $k"

            s = read_struct_from_hdf5(f, T, group=k, use_measurements=use_measurements)
            push!(structs, k => s)
        end
    end

    return structs
end


"""
    read_ordered_structs(filename, T; kwargs...)
"""
function read_ordered_structs(filename::String, T; kwargs...)
    structs = read_structs_from_hdf5(filename, T; kwargs...)
    return structs_to_int_pairs(structs)
end


"""
    structs_to_int_pairs(structs)

Converts a vector of of pairs of string=>struct to a sorted vector of pairs
of int=>struct.
"""
function structs_to_int_pairs(structs::AbstractVector{<:Pair{String, <:Any}})
    structs_new = Pair{Int, <:Any}[(parse(Int, k) => v) for (k, v) in structs]
    idx = sortperm(first.(structs_new))
    structs_new = structs_new[idx]
    return structs_new
end


"""
    write_struct_to_hdf5(filename, obj[, group=""])

Writes a struct to an HDF5 file. Each field in the struct is stored
in the hdf5 file (or group) named after the field in the struct.
"""
function write_struct_to_hdf5(filename::String, obj; kwargs...)
    if isfile(filename)
        rm(filename, force=true)
    end

    h5open(filename, "w") do f
        write_struct_to_hdf5(f, obj; kwargs...)
    end
end



function write_struct_to_hdf5(h5::HDF5.File, obj; group="")
    if group != "" && group ∉ keys(h5)
        HDF5.create_group(h5, group)
    end

    for field in fieldnames(typeof(obj))
        val = getfield(obj, field)
        if eltype(val) <: Measurement
            set_vector!(h5, group * "/" * String(field), middle.(val))
            set_vector!(h5, group * "/" * String(field) * "_em", lower_error.(val))
            set_vector!(h5, group * "/" * String(field) * "_ep", upper_error.(val))
        elseif val isa AbstractVector
            set_vector!(h5, group * "/" * String(field), val)
        elseif field == :annotations
            set_attrs!(h5, val)
        else
            set_vector!(h5, group * "/" * String(field), val)
        end
    end
end


"""
    read_struct_from_hdf5(filename, T[, group="/"])

Reads a struct from an HDF5 file. Expects the struct to be stored
in the group `group` and each field in the hdf5 file (or group)
to be named after the field in the struct.
The type is then called by `T(; kwargs...)`.
"""
function read_struct_from_hdf5(filename::String, T; group="/", use_measurements=:auto)
    local st
    h5open(filename, "r") do f
        st = read_struct_from_hdf5(f, T, group=group, use_measurements=use_measurements)
    end

    return st
end


function read_struct_from_hdf5(h5::HDF5.File, T; group="/", use_measurements=:auto)
    kwargs = Dict{Symbol, Any}()
    for k in keys(h5[group])
        if Symbol(k) ∉ fieldnames(T)
            if use_measurements == false
                @warn "Field $k not found in struct $T"
            elseif endswith(k, "_em") && k[1:end-3] ∈ keys(h5[group])
                use_measurements = true
            elseif endswith(k, "_err") && k[1:end-3] ∈ keys(h5[group])
                use_measurements = true
            else
                # pass
            end
        end
        field = Symbol(k)
        kwargs[field] = get_vector(h5, "$group/$k")
    end

    if use_measurements == :auto
        use_measurements = false
    elseif !(use_measurements isa Bool)
        @error "invalid use_measurements"
    end

    header = get_attrs(h5[group])
    if length(header) > 0
        kwargs[:annotations] = header
    end

    if use_measurements
        kwargs = collapse_errors(kwargs)
    end
    @debug "final kwargs for struct: $(keys(kwargs))"

    return T(; kwargs...)
end

import FITSIO: FITS
import DataFrames: DataFrame
import HDF5

h5open = HDF5.h5open


"""
    read_fits(filename; hdu=2)

Load a FITS file and return a DataFrame using the specified HDU.
"""
function read_fits(filename::String; hdu=2)
    local df
    FITS(filename, "r") do f
        df = DataFrame(f[hdu])
    end

    return df
end


"""
    write_fits(filename, dataframe; overwrite=false, verbose=false)

Write a DataFrame to a FITS file.
"""
function write_fits(filename::String, frame::DataFrame;
        overwrite=false, verbose=false
    )

    if overwrite
        rm(filename, force=true)
    end
    df = to_dict(frame)

    try 
        column_names = ascii.(keys(df))
    catch e
        if isa(e, ArgumentError)
            throw(ArgumentError("Column names must be ASCII"))
        else
            rethrow(e)
        end
    end


    FITS(filename, "w") do f
        write(f, df)
    end

    if verbose
        println("written to $filename")
    end
end



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
    df = Dict(String(name) => frame[:, name] for name in names(frame))

    return df
end



"""
    set_header!(h5f, header)

Sets the header of an HDF5 file using the given dictionary.
"""
function set_header!(h5f::HDF5.File, header::Dict{String,Any})
    if "Header" ∉ keys(h5f)
        HDF5.create_group(h5f, "Header")
    end
    h5_header = h5f["Header"]
    for (key, val) in header
        set_header_attr(h5f, key, val)
    end
end


"""
    set_header_attr(h5f, key, val)

Sets an attribute in the header of an HDF5 file.
"""
function set_header_attr(h5f::HDF5.File, key::String, val)
    header = HDF5.attrs(h5f["Header"])
    header[key] = val
end



"""
    get_header(h5f)

Returns the header of an HDF5 file as a dictionary.
"""
function get_header(h5f::HDF5.H5DataStore)
    return Dict(HDF5.attrs(h5f["Header"]))
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
function write_structs_to_hdf5(filename::String, structs::Vector, labels=nothing)
    if isfile(filename)
        rm(filename, force=true)
    end

    if labels === nothing
        labels = [string(i) for i in 1:length(structs)]
    end

    h5open(filename, "w") do f
        for (label, st) in zip(labels, structs)
            HDF5.create_group(f, label)
            for field in fieldnames(typeof(st))
                val = getfield(st, field)
                set_vector!(f, label * "/" * String(field), val)
            end
        end
    end

end



"""
    read_structs_from_hdf5(filename, T)

Reads a vector of structs from an HDF5 file.
""" 
function read_structs_from_hdf5(filename::String, T)
    structs = T[]
    labels = String[]

    h5open(filename, "r") do f
        for k in keys(f)
            kwargs = Dict{Symbol, Any}()
            for field in fieldnames(T)
                kwargs[field] = get_vector(f, k * "/" * String(field))
            end
            push!(labels, k)

            push!(structs, T(; kwargs...))
        end
    end

    return structs, labels
end

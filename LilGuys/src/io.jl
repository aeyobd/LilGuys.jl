import FITSIO: FITS
import DataFrames: DataFrame
import HDF5

h5open = HDF5.h5open


"""
    load_fits(filename; hdu=2)

Load a FITS file and return a DataFrame using the specified HDU.
"""
function load_fits(filename::String; hdu=2)
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

    FITS(filename, "w") do f
        write(f, df)
    end

    if verbose
        println("written to $filename")
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



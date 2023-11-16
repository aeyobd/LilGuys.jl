using HDF5


# Make a default header for an HDF5 file
function make_default_header(N, mass)
    header = Dict{String,Any}()
    header["NumPart_ThisFile"] = [0, N]
    header["NumPart_Total"] = [0, N]
    header["MassTable"] = [0.0, mass]
    header["Time"] = 0.0
    header["Redshift"] = 0.0
    header["BoxSize"] = 350.0
    header["NumFilesPerSnapshot"] = 1
    return header
end

# Set the header in an HDF5 file
function set_header!(h5f::HDF5.File, header::Dict{String,Any})
    if "Header" ∉ keys(h5f)
        create_group(h5f, "Header")
    end
    h5_header = h5f["Header"]
    for (key, val) in header
        set_header_attr(h5f, key, val)
    end
end


# gets the gadget header of an HDF5 file
function get_header(h5f::HDF5.File)
    return Dict(attrs(h5f["Header"]))
end

# gets a vector from an HDF5 file
function get_vector(h5f::HDF5.File, key::String; mmap=false, group="PartType1")
    path = group * "/" * key
    if mmap
        return HDF5.readmmap(h5f[path])
    else
        return read(h5f[path])
    end
end

function set_vector!(h5f::HDF5.File, key::String, val, group="PartType1")
    if group ∉ keys(h5f)
        create_group(h5f, group)
    end
    path = group * "/" * key
    h5f[path] = val
end


function set_vector_ele!(h5f::HDF5.File, key::String, el::Int, val, group="PartType1")
    path = group * "/" * key
    h5f[path][el] = val
end

function set_header_attr(h5f::HDF5.File, key::String, val)
    header = attrs(h5f["Header"])
    header[key] = val
end



module HDF5Utils
export make_default_header, set_header!, get_header
export set_vector!, get_vector
export get_epsilon


using HDF5
using ..Units


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
function get_vector(h5f::HDF5.File, key::String)
    return read(h5f["PartType1/" * key])
end

function set_vector!(h5f::HDF5.File, key::String, val)
    if "PartType1" ∉ keys(h5f)
        create_group(h5f, "PartType1")
    end
    h5f["PartType1/" * key] = val
end

function set_header_attr(h5f::HDF5.File, key::String, val)
    header = attrs(h5f["Header"])
    header[key] = val
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


end # module HDF5Utils

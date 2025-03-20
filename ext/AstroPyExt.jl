module AstroPyExt
using OrderedCollections

import DataFrames: DataFrame

import LilGuys: read_fits, write_fits

using PythonCall
const AstroPyTable = Ref{Py}()
const Numpy = Ref{Py}()

function __init__()
    AstroPyTable[] = pyimport("astropy.table")
    Numpy[] = pyimport("numpy")
end


"""
    read_fits(filename; hdu=2)

Load a FITS file and return a DataFrame using the specified HDU.

NOTE. This function used to use FITSIO. However, this package 
still is not fully mature and has a tendency to segfault do to poor management
of c-pointers (not easy in julia).
Now, fits are read in using astropy.
"""
function read_fits(filename::String; hdu=2, columns=nothing)
    df = DataFrame()
    np = Numpy[]

    table = AstroPyTable[].Table.read(filename; hdu=hdu, format="fits")

    all_columns = pyconvert(Vector{String}, table.columns)

    if isnothing(columns)
        columns = all_columns
    else
        if !(columns ⊆ all_columns)
            @error "Columns not found in data: $(setdiff(columns, all_columns))"
        end
    end

    for colname in pyconvert(Vector{String}, (columns))
        coldat = np.array(table[colname])
        @debug "reading $colname"
        @debug "col[0] $(coldat[0])"
        dtype = coldat.dtype
        local x

        if pyconvert(Bool, np.issubdtype(dtype, np.bytes_))
            coldat = coldat.astype(np.str_)
            @debug "reading as bytes"
            x = pyconvert(Vector{String}, coldat)
        elseif pyconvert(Bool, np.issubdtype(dtype, np.floating))
            @debug "reading as float"
            x = pyconvert(Vector{Float64}, coldat)
        elseif pyconvert(Bool, np.issubdtype(dtype, np.integer))
            @debug "reading as float"
            x = pyconvert(Vector{Int}, coldat)
        else
            @debug "reading as other"
            x = pyconvert(Vector, coldat)
        end

        df[!, colname] = x
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

    df = OrderedDict(
        col => Numpy[].array(frame[:, col]) for col in names(frame)
       )

    for col in names(frame)
        if eltype(frame[!, col]) <: AbstractString
            df[col] = Numpy[].array(frame[:, col], Numpy[].str_)
        end
    end


    try 
        column_names = ascii.(keys(df))
    catch e
        if isa(e, ArgumentError)
            throw(ArgumentError("Column names must be ASCII"))
        else
            rethrow(e)
        end
    end


    table = AstroPyTable[].Table(PyDict(df))
    table.write(filename)

    if verbose
        println("written to $filename")
    end
end


end # module

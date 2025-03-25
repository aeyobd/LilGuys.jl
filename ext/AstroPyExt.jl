module AstroPyExt
using OrderedCollections

import DataFrames: DataFrame

import LilGuys: read_fits, write_fits

using PythonCall
const AstroPyTable = Ref{Py}()
const Numpy = Ref{Py}()
const Pandas = Ref{Py}()


function __init__()
    AstroPyTable[] = pyimport("astropy.table")
    Numpy[] = pyimport("numpy")
    Pandas[] = pyimport("pandas")
end


function pycol_to_vec(col)
    dtype = col.dtype
end

"""
    read_fits(filename; hdu=2)

Load a FITS file and return a DataFrame using the specified HDU (1-indexed).

NOTE. This function used to use FITSIO. However, this package 
still is not fully mature and has a tendency to segfault do to poor management
of c-pointers (not easy in julia).
Now, fits are read in using astropy.
"""
function read_fits(filename::String; hdu=2, columns=nothing, kwargs...)
    df = DataFrame()

    tab = AstroPyTable[].Table.read(filename; hdu=hdu-1, format="fits", kwargs...)
    table = tab.to_pandas()
    @info "pandas table loaded in"

    all_columns = pyconvert(Vector{String}, table.columns)

    if isnothing(columns)
        columns = all_columns
    else
        if !(columns ⊆ all_columns)
            @error "Columns not found in data: $(setdiff(columns, all_columns))"
        end
    end

    size_filter = [length(table[col].shape) <= 1 for col in columns]
    if sum(.!size_filter) > 0
        @warn "skipping multidimensional columns: $(columns[.!size_filter])"
    end
    columns = columns[size_filter]



    for colname in columns
        @debug colname
        @debug table[colname].dtype
        @debug table[colname].shape

        df[!, colname] = read_pandas_column(table[colname])
        #table.drop(columns = colname) # free memory as we go
    end

    return df
end


"""
    read_pandas_column(column)

Reads in a pandas series 
"""
function read_pandas_column(column)
    if pyisinstance(column.dtype, Numpy[].dtypes.ObjectDType) && hasproperty(column, :str)
        @debug "reading as bytes"
        col = column.str.decode("ascii") # fits is ascii only
        coldat = map(col) do x
            if pyisinstance(x, pybuiltins.str)
                return pyconvert(String, x)
            elseif pyconvert(Bool, pybool(Numpy[].isnan(x)))
                return NaN
            else
                @error "type $(pybuiltins.type(x)) not implemented for bytes"
            end
        end

    else
        coldat = pyconvert(Vector, column)
    end

    return coldat
end


"""
    write_fits(filename, dataframe; overwrite=false, verbose=false)

Write a DataFrame to a FITS file.
"""
function write_fits(filename::String, frame::DataFrame;
        overwrite=false, verbose=false, kwargs...
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


    table = AstroPyTable[].Table(PyDict(df); kwargs...)

    for col in table.columns.keys()
        @debug col
        @debug table[key].dtype
        @debug table[key].shape
    end

    table.write(filename)

    if verbose
        println("written to $filename")
    end
end


end # module

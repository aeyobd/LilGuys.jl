import FITSIO: FITS
import DataFrames: DataFrame


function load_fits(filename::String; hdu=2)
    local df
    FITS(filename, "r") do f
        df = DataFrame(f[hdu])
    end

    return df
end


function write_fits(filename::String, frame::DataFrame;
        overwrite=false, verbose=false
    )

    if overwrite
        rm(filename, force=true)
    end
    df = Dict(frame)

    FITS(filename, "w") do f
        write(f, df)
    end

    if verbose
        println("written to $filename")
    end
end


function Dict(frame::DataFrame)
    df = Dict(String(name) => frame[:, name] for name in names(frame))

    return df
end

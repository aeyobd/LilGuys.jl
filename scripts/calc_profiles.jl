import LilGuys as lguys
using ArgParse
using CSV
using DataFrames


function main(dirname, verbose=true)
    out = lguys.Output(dirname)
    cens = lguys.cens(out)
    for i in 1:length(cens)
        prof = lguys.Profile(out[i], cens[i].x_c, cens[i].v_c)
        lguys.write("profiles/$(i).csv", prof)

        if verbose
            println("Wrote profiles/$(i).csv")
        end
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    ap = ArgParse.ArgumentParser()
    ArgParse.add_arg!(ap, "dirname"; help="Directory name")
    args = ArgParse.parse(ap)
    println("reading from $(args.dirname)")
    main(args.dirname)
end


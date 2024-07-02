struct AgamaPotential
    _ptr::Ptr{CVoid}
end


function AgamaPotential(;kwargs...) 
    if length(kwargs) < 1
        ArgumentError("Agama requires at least one kwarg")
    end


    params = string(kwargs)
    params = params[2:end-1] # strip parenthasis

    ptr = @ccall "./agama.so".agama_createPotential(args::Cstring)::Ptr{Cvoid}

    return AgamaPotential(ptr)
end




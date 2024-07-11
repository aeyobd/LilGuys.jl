
struct AgamaPotential
    _ptr::Ptr{Cvoid}
end

agama_path = "agama"

function AgamaPotential(; kwargs...) 
    if length(kwargs) < 1
        ArgumentError("Agama requires at least one kwarg")
    end


    args = string(NamedTuple(kwargs))
    args = args[2:end-1] # strip parenthasis
    args = replace(args, "\""=>"") # no quotes
    args = replace(args, " "=>"") # no quotes
    if length(kwargs) == 1
        args = args[1:end-1] # strip trailing comma
    end


    println("args: ", args)
    ptr = @ccall agama_path.agama_createPotential(args::Cstring)::Ptr{Cvoid}


    if ptr == C_NULL
        raise_agama_error()
    end

    return AgamaPotential(ptr)
end

function raise_agama_error()
    message = @ccall agama_path.agama_getError()::Cstring
    message = unsafe_string(message)
    error("Agama error: $message")
end



function get_Φ(pot::AgamaPotential, pos::Vector{Float64}, time::Float64)
    deriv = zeros(3)
    deriv2 = zeros(3)

    result = @ccall agama_path.agama_evalPotential(pot._ptr::Ptr{Cvoid}, pos::Ptr{Cdouble}, time::Cdouble, deriv::Ptr{Cdouble}, deriv2::Ptr{Cdouble})::Cdouble

    return result
end

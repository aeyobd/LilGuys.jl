"""The path to agama c library"""
_agama = "agama"


"""
A wrapper around a Agama potential C object.
"""
struct AgamaPotential
    _ptr::Ptr{Cvoid}
end



function AgamaPotential(; kwargs...) 
    args = _kwargs_to_Cstring(; kwargs...)

    ptr = @ccall _agama.agama_createPotential(args::Cstring)::Ptr{Cvoid}

    if ptr == C_NULL
        raise_agama_error()
    end

    # need preserve to avoid segfault
    return @GC.preserve AgamaPotential(ptr)
end



function raise_agama_error()
    message = @ccall _agama.agama_getError()::Cstring
    message = unsafe_string(message)
    error("Agama error: $message")
end


"""
    get_Φ(pot, pos, time)

Compute the potential at a given position and time.
"""
function get_Φ(pot::AgamaPotential, pos::Vector{Float64}, time::Float64)
    deriv = zeros(3)
    deriv2 = zeros(3)

    result = @ccall _agama.agama_evalPotential(pot._ptr::Ptr{Cvoid}, pos::Ptr{Cdouble}, time::Cdouble, deriv::Ptr{Cdouble}, deriv2::Ptr{Cdouble})::Cdouble

    return result
end


"""
    get_ρ(pot, pos[, time])

Compute the density at a given position and time.
"""
function get_ρ(pot::AgamaPotential, pos::Vector{Float64}, time::Float64=0.)
    result = @ccall _agama.agama_evalDensity(pot._ptr::Ptr{Cvoid}, pos::Ptr{Cdouble}, time::Cdouble)::Cdouble

    return result
end

function get_R_circ(pot::AgamaPotential, E::Float64)
    result = @ccall _agama.agama_R_circ(pot._ptr::Ptr{Cvoid}, E::Cdouble)::Cdouble

    return result
end


function R_from_Lz(pot::AgamaPotential, Lz::Float64)
    result = @ccall _agama.agama_R_from_Lz(pot._ptr::Ptr{Cvoid}, Lz::Cdouble)::Cdouble

    return result
end


function get_R_max(pot::AgamaPotential, E::Float64)
    result = @ccall _agama.agama_Rmax(pot._ptr::Ptr{Cvoid})::Cdouble

    return result
end



function _kwargs_to_Cstring(; kwargs...)
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
    return args
end


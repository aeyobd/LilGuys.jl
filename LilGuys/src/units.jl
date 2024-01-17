
export G, M0, R0, T0, V0, F

const G = 1.
const M0 = 1e10 #Msun
const R0 = 1. # kpc
const T0 = 4.718e6 # years
const V0 = 207.4 # km/s

const G_U = PyNULL()
const M_U = PyNULL()
const R_U = PyNULL()
const T_U = PyNULL()
const V_U = PyNULL()

const G_AP = PyNULL()
const M_AP = PyNULL()
const R_AP = PyNULL()
const T_AP = PyNULL()
const V_AP = PyNULL()

F = Float64
OptVector = Union{Vector{F}, Nothing}
OptMatrix = Union{Matrix{F}, Nothing}

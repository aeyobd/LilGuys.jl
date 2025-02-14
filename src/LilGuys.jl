module LilGuys

export Snapshot
export Output
export ICRS, PhasePoint
export to_galcen, to_sky
export save

export calc_r, calc_v
export @savefig


# profile tools
export AbstractProfile, NFW
export calc_ρ, calc_M, calc_r_circ_max, calc_v_circ_max, calc_v_circ

export M2MSUN, R2KPC, V2KMS, T2GYR
export ICRS, HelioRest, Galactocentric, transform



include("units.jl")
include("interface.jl")
using .Interface


include("utils.jl")

include("coordinates.jl")

include("io.jl")
include("snapshot.jl")
include("output.jl")

include("spherical.jl")
include("coord_trans.jl")   

include("physics.jl")
include("gravity.jl")

include("analytic_profiles.jl")
include("nfw.jl")
include("scaling_relations.jl")

include("project.jl")
include("density_3d.jl")
include("density_2d.jl")
include("density_3d_star.jl")

include("orbits.jl")
include("potentials.jl")

include("centres/Centres.jl")


# Empty function definitions for MakieExt
public plot_xyz, plot_xyz!
public projecteddensity, projecteddensity!
public hide_grid!, cmd_axis

function plot_xyz end
function plot_xyz! end
function cmd_axis end

function projecteddensity end
function projecteddensity! end
function hide_grid! end

macro savefig end



using .Centres

end # module LilGuys

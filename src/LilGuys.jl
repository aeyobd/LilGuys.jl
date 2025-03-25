module LilGuys



export M2MSUN, R2KPC, V2KMS, T2GYR
public F, G
public arcmin2kpc, kpc2arcmin, pm2kms, kms2pm, dm2dist, dist2dm
include("units.jl")

export Measurement, middle, lower_bound, upper_bound, credible_interval
include("measurements.jl")


public mean, std, midpoints, quantile, variance
public expinti, erf
public integrate, histogram 
public effective_sample_size, bins_default, bins_equal_number, bins_both, bins_equal_width, default_bin_width, default_n_per_bin
public DataFrame

include("interface.jl")
using .Interface

public ConstVector
public struct_to_dict, dict_to_tuple
public collapse_errors
public randu, rand_unit
public gradient
public lerp
public normal_cdf, gaussian, logistic, logit
public centroid, centroid_err
include("utils.jl")

public CoordinateFrame, AbstractCartesian, AbstractSkyCoord
public Point3D, Point6D, SkyCoord
public GalactocentricFrame
export ICRS, GSR, Galactocentric, Cartesian
public position, velocity
public rand_coord, rand_coords, coord_from_file, coords_from_df
include("coordinates.jl")

export read_hdf5_table, write_hdf5_table
export read_struct_from_hdf5, read_structs_from_hdf5
export write_struct_to_hdf5, write_structs_to_hdf5
include("io.jl")

public write
export Snapshot
include("snapshot.jl")


export Output
public extract, extract_vector, peris_apos
include("output.jl")


public to_tangent, angular_distance, unit_vector
public cartesian_to_sky, rotate_sky, Rx_mat, Ry_mat, Rz_mat
include("spherical.jl")
public transform
include("coord_trans.jl")   


export radii, speed
public L_tot, L_spec, calc_E_tot, calc_E_spec
public calc_ϵ, calc_K_spec, calc_K_tot, calc_W_tot
public x_velocity, y_velocity, z_velocity, x_position, y_position, z_position
public bound_particles, bound_particles_recursive_1D
include("physics.jl")

public DistributionFunction
public potential_spherical, potential_spherical_func, potential_nbody, potential_spherical_discrete
public F_grav
include("gravity.jl")

# profile tools
export AbstractProfile, NFW
public Plummer, KingProfile, Exp2D, Exp3D
export density, surface_density, mass, v_circ
include("analytic_profiles.jl")

export r_circ_max, v_circ_max
include("nfw.jl")

include("scaling_relations.jl")

include("project.jl")

include("density_2d.jl")

export MassProfile3D
include("density_3d.jl")

include("density_3d_star.jl")

# orbits
export Orbit
public to_frame, resample, leap_frog
include("orbits.jl")

include("potentials.jl")

include("centres/Centres.jl")
using .Centres



# Empty function definitions for MakieExt

public plot_xyz, plot_xyz!
public plot_log_Σ!, plot_Γ!
public projecteddensity, projecteddensity!
public hide_grid!, cmd_axis
export @savefig


function plot_xyz end
function plot_xyz! end
function plot_log_Σ! end
function plot_Γ! end
function cmd_axis end

function projecteddensity end
function projecteddensity! end
function hide_grid! end

macro savefig end


# Empty functions for AgamaExt
public AgamaPotential
public sample_potential
function AgamaPotential end
function sample_potential end


# empty functions for AstroPyExt
public read_fits, write_fits
"""
import PythonCall to use this method
"""
function read_fits end

"""
import PythonCall to use this method
"""
function write_fits end


end # module LilGuys

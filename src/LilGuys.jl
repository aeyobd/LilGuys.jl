module LilGuys

export 
    M2MSUN, 
    R2KPC, 
    V2KMS, 
    T2GYR

public 
    F, G,
    arcmin2kpc, kpc2arcmin, 
    pm2kms, kms2pm, 
    dm2dist, dist2dm

include("units.jl")


export 
    Measurement, 
    middle, 
    lower_error, 
    upper_error, 
    error_interval

include("measurements.jl")


public 
    DataFrame,

    mean, 
    std, 
    quantile, 
    variance,
    midpoints, 

    expinti, 
    erf,

    integrate, 
    curve_fit,

    histogram, 
    bins_default,
    bins_equal_number,
    bins_both,
    bins_equal_width,
    default_bin_width,
    default_n_per_bin


include("interface.jl")
using .Interface


public 
    ConstVector,
    struct_to_dict, 
    dict_to_tuple,
    collapse_errors,
    randu, rand_unit,
    gradient,
    lerp,
    normal_cdf, 
    gaussian, 
    logistic, 
    logit,
    centroid, 
    centroid_err,
    effective_sample_size

include("utils.jl")


public 
    CoordinateFrame, 
    AbstractCartesian, 
    AbstractSkyCoord,
    Point3D, 
    Point6D, 
    SkyCoord,
    GalactocentricFrame,
    position, 
    velocity, 
    rand_coord,
    rand_coords,
    coord_from_file,
    coords_from_df

export 
    ICRS, 
    GSR, 
    Galactocentric, 
    Cartesian

include("coordinates.jl")


export 
    read_hdf5_table, 
    write_hdf5_table,
    read_struct_from_hdf5, 
    read_structs_from_hdf5,
    write_struct_to_hdf5, 
    write_structs_to_hdf5

include("io.jl")


public write
export Snapshot

include("snapshot.jl")


export Output
public 
    extract, 
    extract_vector, 
    peris_apos

include("output.jl")


public 
    to_tangent, 
    angular_distance, 
    unit_vector,
    cartesian_to_sky, 
    rotate_sky, 
    Rx_mat, 
    Ry_mat, 
    Rz_mat

include("spherical.jl")


public transform

include("coord_trans.jl")   


export radii, speed
public 
    L_tot, 
    L_spec, 
    calc_E_tot, 
    calc_E_spec,
    calc_ϵ, 
    calc_K_spec,
    calc_K_tot,
    calc_W_tot,
    x_velocity,
    y_velocity,
    z_velocity,
    x_position,
    y_position,
    z_position,
    bound_particles,
    bound_particles_recursive_1D

include("physics.jl")


public 
    DistributionFunction,
    potential_spherical,
    potential_spherical_func,
    potential_nbody,
    potential_spherical_discrete,
    F_grav

include("gravity.jl")


# profile tools
export 
    AbstractProfile, 
    NFW,
    density, 
    surface_density, 
    mass, 
    v_circ

public 
    Plummer, 
    KingProfile, 
    Exp2D, 
    Exp3D

include("analytic_profiles.jl")

export r_circ_max, v_circ_max
public R200, M200

include("nfw.jl")


include("scaling_relations.jl")


public 
    calc_R_ell, 
    calc_R_ell_sky, 
    shear_points_to_ellipse,
    calc_centre2D, 
    to_orbit_coords,
    to_gaia

include("project.jl")


export StellarDensityProfile
public filter_empty_bins

include("density_2d.jl")


export MassProfile3D

include("mass_profile_3d.jl")

export DensityProfile3D
public radial_velocities

include("density_3d.jl")


export StellarDensity3D

include("stellar_density_3d.jl")


# orbits
export Orbit
public to_frame, resample, leap_frog

include("orbits.jl")

include("potentials.jl")

include("centres/Centres.jl")
using .Centres



# Empty function definitions for MakieExt

export @savefig
public 
    plot_xyz,
    plot_xyz!,
    plot_log_Σ!, 
    plot_Γ!, 
    plot_log_Σ, 
    plot_Γ,
    projecteddensity, 
    projecteddensity!,
    hide_grid!, 
    cmd_axis



function plot_xyz end
function plot_xyz! end
function plot_log_Σ! end
function plot_log_Σ end
function plot_Γ! end
function plot_Γ end
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



end # module LilGuys

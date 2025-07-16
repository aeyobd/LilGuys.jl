# LilGuys
Dwarf galaxies are just 'lilguys. 

## Notes about structure and conventions

### Units
All internal code is in special code units, where G = 1, M = 1, R = 1, and t = 1.
the units file will contain methods to covert to physical units

### Objects

1. `Point` - a point in 3D space
2. `PhasePoint` - a point in 6D phase space
3. `Particle` - points
4. `ICRS`, `GSR`, `Cartesian`, `Galactocentric`: Sky coordinate types
5. `StellarDensityProfile`, `StellarMassProfile`
6. `DensityProfile`, `StellarDensityProfile3D`, `MassProfile`


### Strucutre
- `units.jl` A few basic unit conversions and our code unit definitions
- `utils.jl`
- `interface.jl`
- `measurements.jl`
- `io.jl`
  - `snapshot.jl`
  - `output.jl`

- `coordinates.jl`: contains methods to convert between galactocentric and geocentric frames
  - `coord_trans.jl`
  - `project.jl`
  - `spherical.jl`

- `gravity.jl`: gravitational methods: profiles, forces, and potentials
- `physics.jl`
- `stellar_density_3d.jl`
- `density_2d.jl`
- `density_3d.jl`
- `analytic_profiles.jl`
- `nfw.jl`
- `scaling_relations.jl`
- `potentials.jl`, `orbits.jl`


## Naming conventions and types

### Variables
Julia is column major. As such, vectors are columns, so if we are working in R^3, then `\mathbf{x}` is represented as a 3-vector, and a set of N points is represented as a 3xN matrix. 
Technically, scalar fields should be 1xN matrices (as julia vectors are Nx1), which is something I would like to move to but haven't completed yet

For single scalars (mathematically $\in \R$), I like single letters like $r$ (radius), $\Phi$ (potential), $s$ (scalar velocity), $m$ (mass) and so on. Especially as this is a mathematical project, I think this is okay for now, but more verbose names can be added later.


### Functions

Common functions (sometimes imported)  like `mean`, `std`, and `norm` are noun clauses.



### Examples

Observables (below) may also specify symmetric errors as `_err` or assymetric errors with `_em` and `_ep`.  However, these are only used in stored files, otherwise represented by internal `Measurement` type

- `ra`: Right ascension (ICRS J2000) in degrees
- `dec`: Declination (ICRS J2000) in degrees
- `distance`: heliocentric distance in kpc
- `pmra`: Corrected proper motion in right ascension (mas/yr). Name is technically inconsistent with underscore rule but kept for consistency with Gaia
- `pmdec`: Corrected proper motion in declination (mas/yr).
- `radial_velocity`: Heliocentric line-of-sight radial velocity of a stars (km/s). 
- `position_angle`: Position angle of satillite major axis in degrees (North to East). 
- `ellipticity`: Ellipticity of satellite density profile (1 - b/a?)
- `xi`, `eta`: tangent plane coordinates (ra / dec) in arcminutes. Based on `ra` and `dec` and presently assumed centre (probably from `observations/galaxyname/observed_properties.toml` or simulation centre)
- `R_h`: Half light radius of satellite in *arcminutes*
- `R_ell`: elliptical radius in units of arcminutes
- `sigma_v`: LOS velocity dispersion (in arcminutes)



Generic quantities

lengths

- `log` as prefix, always log 10
- `R` 2D radius / projected radius (vector)
- `radius`
- `radii` (`r` internally) 3D radius vector
- `log_R`, `log_radii` like above but log10
- `position`, `positions` 3xN vector of positions
- `r_h`
- `R_h`
- `break_radius`

speeds

- `velocities` 3xN vector of velocities 
- `v_circ` circular velocity (from enclosed mass)

acceleration / force

- `accelerations`
- `force`

masses

- `mass` total or enclosed mass
- `masses` point particle masses
- `mass_enclosed`
- `M` (attribute/scale)

2d densities

- `surface_density`
- `log_surface_density`

3D densities

- `density`
- `log_density`

energy

- `energy`
- `binding_energy`
- `surface_density`
- `potential`, `potential_ext` Gravitational potential energy $\Phi$.



momentum

- `angular_momentum`



# Development

![image-20250512164547349](/Users/daniel/Library/Application Support/typora-user-images/image-20250512164547349.png)

Clean and basic tests

- units (100%)
- utils (93%) 
- interface (94%)

TO check

- io 
- snapshot
- output
- coord_trans
- spherical
- physics
- gravity
- analytic_profiles
- nfw
- scaling_relations
- project
- Centres.jl
- density_2d
- measurements (basic)
- MassQuantiles
- MassWithinRadii

Needs tests:

- mass_profile_3d
- MakieExt
- *Shrinking spheres robust tests*
- Coordinate constructors more tests
- cylindrical coordinate frame & tests
- velocity anisotropy first tests

More features

- Create working density centres (fuzzy centres)
- First tests for potentials
- first tests for orbits
- more tests for project.jl
- FITS metadata

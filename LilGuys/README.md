# LilGuys
LIttLe Galaxies Undergoing milkY way Stripping

## Notes about structure and conventions

### Units
All internal code is in special code units, where G = 1, M = 1, R = 1, and t = 1.
the units file will contain methods to covert to physical units

### Objects

1. `Point` - a point in 3D space
2. `PhasePoint` - a point in 6D phase space
2. `Particle` - points


### Strucutre
- units.j
- coordinates.jl: contains methods to convert between galactocentric and geocentric frames
- gravity.jl: gravitational methods: profiles, forces, and potentials
- 


## Naming conventions and types

### Variables
Julia is column major. As such, vectors are columns, so if we are working in R^3, then `\mathbf{x}` is represented as a 3-vector, and a set of N points is represented as a 3xN matrix. 
Technically, scalar fields should be 1xN matrices (as julia vectors are Nx1), which is something I would like to move to but haven't completed yet

For single scalars (mathematically $\in \R$), I like single letters like $r$ (radius), $\Phi$ (potential), $s$ (scalar velocity), $m$ (mass) and so on. Especially as this is a mathematical project, I think this is okay for now, but more verbose names can be added later.

As the two most common vectors ($\mathbf{x}$ and $\mathbf{v}$) have overloaded names (V is potential, x is a coordinate), I will use three letter abbreviations for 
- `x_vec`
- `v_vec`
- `a_vec`
pluralizing to `x_vecs`, for 3xN matrices of points, etc.
Subscripts are added after to (`x_vec_i` pluralizing to `x_vecs_i` so dimensionality are contained in the same group (vecs)).


### Functions

Function names are snake case, and for the most part should be verb clauses. This helps avoid accidentally overriding functions which represent common variables like radius `r`. So, most physical quantity calculations begin with `calc_r`, etc. 

Common functions (sometimes imported)  like `mean`, `std`, and `norm` are noun clauses but these are exceptions in this project.

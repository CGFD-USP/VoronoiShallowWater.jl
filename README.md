# VoronoiShallowWater.jl

`VoronoiShallowWater.jl` is a Shallow Water model that uses Voronoi Tesselation meshes.
Spherical and planar bi-periodic meshes are supported.
For planar bi-periodic Voronoi meshes creation see [VoronoiMeshes.jl](https://github.com/favba/VoronoiMeshes.jl)

## Install Guide 

Assuming Julia already installed (for instance via juliaup).

Using Julia version 1.11 or higher:
```julia
import Pkg
Pkg.activate(temp=true)
Pkg.add(url="https://github.com/CGFD-USP/VoronoiShallowWater.jl.git")

# Optional packages, for grid creation; plotting; import / export to NetCDF; and import / export to VTK.
Pkg.add(url="https://github.com/favba/VoronoiMeshes.jl.git")
Pkg.add("DelaunayTriangulation")
Pkg.add("GLMakie")
Pkg.add("NCDatasets")
Pkg.add("ReadVTK")
Pkg.add("WriteVTK")
```

If using Julia v1.10, then the unregistered dependencies must be explicitly installed before installing this package:
```julia
import Pkg
Pkg.activate(temp=true)
Pkg.add(url="https://github.com/favba/TensorsLite.jl.git")
Pkg.add(url="https://github.com/favba/TensorsLiteGeometry.jl.git")
Pkg.add(url="https://github.com/favba/VoronoiMeshes.jl.git")
Pkg.add(url="https://github.com/favba/VoronoiOperators.jl.git")
Pkg.add(url="https://github.com/CGFD-USP/VoronoiShallowWater.jl.git")

# Optional packages, for grid creation; plotting; import / export to NetCDF; and import / export to VTK.
Pkg.add("DelaunayTriangulation")
Pkg.add("GLMakie")
Pkg.add("NCDatasets")
Pkg.add("ReadVTK")
Pkg.add("WriteVTK")
```
The package and the dependencies will be ported into registered packages soon.

## Quick usage guide

The simulation is configured by a single TOML file and the simulation is run, in a julia section with:
```julia
julia> using VoronoiShallowWater

julia> run_simulation("config-example.toml")
```

Or, from the command line, for an optimized build:
```bash
julia -O3 --threads auto --project=@where_instaled -e 'using VoronoiShallowWater; run_simulation("config-example.toml")'
```

A sample TOML configuration file is given below:
```toml
# config-example.toml

model = "TRiSK"  # "Peixoto" is the other option. All String values case insensite in the configuration.
grid_file = "mesh_periodic_regular_nc256.vtu"  # Or NetCDF file "mesh_name.nc", any file readable by `VoronoiMesh`.
grid_scale_factor = 4.003155589280872e7  # This is 2*pi*radius_earth. Default is not to perform any scaling.
time_step = { CFL = 0.14, maximum_velocity = 100.0 }  # dt = CFL * dx/(sqrt(H*g) + maximum_velocity) or float for time in seconds.
duration = "20_00:00:00"  # Or float (time in seconds)

# initial_condition = "Zonal geostrophic balance" # Or table like below for changing default values or string with netcdf file name.
[initial_condition]
name = "Zonal geostrophic balance"
H = 1e4  # Mean Height
u0 = 50.0  # Zonal jet maximum velocity
g = 9.80616  # Gravity acceleration
f = 1.4584e-4  # Constant Coriolis parameter

[output]
file_name = "output.nc" # NetCDF file name where the save the results
time = "01_00:00:00"  # Save every `time` time. Or float for time in seconds
variables = ["h", "u", "velocity", "pv", "ke"]  # Vector of Variables to save to file. Default is ["h", "u"]
precision = "Float32"  # Same as "single". Save variables in `precision` precision. Default is "Float64" ("double"). Also accepts "Float16" ("half")

# Print a summary of the simulation during execution.
[summary]
enabled = true
steps = 20  # Write summary every `steps` steps.
min_max = true  # print min and max of `u` and `h`.
energy_terms = true  # print the mean kinetic and available potential energy terms

```


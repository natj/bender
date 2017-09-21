# Bender


### Short code dissection and structure

* `main`: Main loop for calling sub-functions in correct order
* `bender`: Hamilton-Jacobi ray tracing propagator
* `beaming`: Construction of the beaming function; includes code for solving the electron scattering atmosphere
* `img`: Computes the image visible for the observer.
* `cspot`: Time-dependent Cubature integration of the flux from the observers image; used for pulse profile calculations
* `radiation`: Handles actual emerging spectra calculations; currently Planck function 
* `rtrace`: Stretched polar image plane pixelation that constructs interpolators for flux and other quantities
* `rtrace_cart`: Cartesian image plane pixelation; *Use `rtrace` instead*.
* `strig`: Auxiliary (spherical) trigonometric functions


##### Additional tools & scripts
* `3d_paths`: Short script for saving full photon paths.
* `datafy`: *very* simple synthetic data constructor script for pulse profiles
* `skymap`: Script for looping over different observers for computing skymaps.
* `src/`: Includes a crude sketch of a C++ version. **TODO**

# Bender
### Radiation from rapidly rotating oblate neutron stars
----------

Bender is a ray tracing code for computing radiation from rotating neutron stars. 
Or to put it in another way: Bender is a collection of scripts written in high-level programming language for easy adaptation and modification.
We take into account all of the important special and general relativistic effects that change the emerging radiation. 
Hence, it can be used to compute the observed energy spectra of emerging radiation when the initial flux energy distribution is known. 
Additionally, we can deal with different beaming patterns, i.e., angle dependency of the radiation from the NS surface. 
There is also a possibility to compute time-dependent spectra for e.g., pulse profile computations. 
For a more in-depth discussion about the physics, see the main paper Nättilä & Pihajoki (2016).

** If you want to use bender for your own calculations, I am happy to help you to get things running!**


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

# Bender
### Radiation from rapidly rotating oblate neutron stars
----------

Bender is a ray tracing code for computing radiation from rotating neutron stars. 
Or to put it in another way: Bender is a collection of scripts written in high-level programming language for easy adaptation and modification.
We take into account all of the important relativistic effects that change the emerging radiation. 
Hence, it can be used to compute the observed energy spectra of emerging radiation when the initial flux energy distribution is known. 
Additionally, we can deal with different beaming patterns, i.e., angle dependency of the radiation from the NS surface. 
There is also a possibility to compute time-dependent spectra for e.g., pulse profile computations. 
For a more in-depth discussion about the physics, see the main paper Nättilä & Pihajoki (2016).

**If you want to use bender for your own calculations, I am happy to help you to get things running!**


## Code Overview

Bender is a collection of scripts only that you can find in this main folder. 
These are mainly done with python that uses the excellent `arcmancer` library that actually propagate the photons.
Hence, bender is actually just a collection of integrators and other such methods to simplify these kind of calculations.
Other option is to use the Hamilton-Jacobi propagator (written in julia) that can be found from julia folder.

First, see:

- `bender.py`: standalone script for computing neutron star image

Then, more complex examples, such as:
- `neutronstar.py`: neutron star image constructed using auxiliary functions and other helper classes
- `pulse.py`: time-dependent pulse profiles from hot spots
- `pulse_polar.py`: same as above but in polar coordinates
- `sweep.py`: sweep over observer inclination angles; useful for computing many angle-dependent static quantities



Auxiliary
- `img.py`: Image plane construction
- `lineprofile.py`: Sum and bin observed image to get line profile
- `radiation.py`: different radiation models
- `spot.py`: spot shapes, sizes, etc.
- `units.py`: transformations from numerical to physical units
- `visualize.py`: visualization scripts 
- `visualize_polar.py`: same as above but in polar grid







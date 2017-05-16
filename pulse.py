import sys
sys.path.append('/Users/natj/projects/arcmancer/lib/')
import pyarcmancer as pyac

from img import Imgplane
from lineprofile import *
from spot import Spot
from visualize import Visualize


import numpy as np
import matplotlib as mpl
from pylab import *

from matplotlib import cm

import scipy.interpolate as interp
from cubature import cubature


#from joblib import Parallel, delayed
#import multiprocessing




mpl.rc('font', family='serif')
mpl.rc('xtick', labelsize='small')
mpl.rc('ytick', labelsize='small')
gs = GridSpec(3, 3)
gs.update(hspace = 0.3)
#gs.update(wspace = 0.3)

#num_cores = multiprocessing.cpu_count()
#print "num of cores {}", num_cores


mpl.rcParams['image.cmap'] = 'inferno'

##################################################
# Star parameters
R = 12.0
M = 1.4
freq = 600.0
incl = 60.0


##################################################
# Spot parameters
rho = 20.0
colat = 45.0


##################################################
# Operate in units where G = c = 1.
# Use units of solar masses
solar_mass_per_km = 0.6772
solar_mass_per_s  = 2.03e5


## conversion factors
kelvin_per_kev = 1.16045e7
km_per_kpc     = 3.086e16
kg_per_Msun    = 1.988435e30

# Variables in units of solar mass are derived here
# and typically presented with full name
mass        = M
radius      = R * solar_mass_per_km / mass
angvel      = freq * 2.0*np.pi / solar_mass_per_s * mass
compactness = np.sqrt(1 - 2/radius)


##################################################
#Create spots
spot = Spot(colat, rho, angvel) #initalize spot(s)



#Setup pyarcmancer
##################################################
conf = pyac.Configuration()

conf.absolute_tolerance       = 1e-12 * radius
conf.relative_tolerance       = 1e-12
conf.henon_tolerance          = 1e-8
conf.sampling_interval        = 1e-3
conf.minimum_stepsize         = 1e-10 * radius
conf.maximum_steps            = 10000
conf.enforce_maximum_stepsize = False
conf.enforce_minimum_stepsize = True
conf.enforce_maximum_steps    = True
conf.store_only_endpoints     = True


#pyac.Log.set_console()
pyac.Log.set_file()

##################################################
#Define metric and surfaces of the spacetime 

metric = pyac.SchwarzschildMetric(mass/radius)
#metric = pyac.AGMMetric(radius, 1.0, angvel, pyac.AGMMetric.MetricType.agm_standard)
#metric = pyac.AGMMetric(radius, 1.0, angvel, pyac.AGMMetric.MetricType.agm_no_quadrupole)

ns_surface = pyac.AGMSurface(radius, 1.0, angvel, pyac.AGMSurface.SurfaceType.spherical)
#ns_surface = pyac.AGMSurface(radius, 1.0, angvel, pyac.AGMSurface.SurfaceType.oblate)
#ns_surface = pyac.AGMSurface(radius, 1.0, angvel, pyac.AGMSurface.SurfaceType.agm)
surfaces = [ ns_surface ]



# Build and configure image plane by hand
img = Imgplane(conf, metric, surfaces)

img.verbose  = 1
img.incl     = np.deg2rad(incl) #set inclination
img.distance = 100.0*mass #set distance


#Locate star edges
img.find_boundaries()

#Build internal coarse grid for the interpolation routines
img.generate_internal_grid(Nrad = 30, Nchi = 30 )
img.dissect_geos()


###################################################
# flux function for integration
def flux(xy):
    x, y = xy
    time, phi, theta, cosa, reds = img.get_pixel(x, y) 
    coords = np.array([time, theta, phi])

    hit = spot.hit(coords)

    if hit:
        return np.array([1.0, 1.0])
    else:
        return np.array([0.0, 0.0])




#Construct output xy image plane from img object
##################################################
ion()
visz = Visualize()
visz.compactness = compactness
visz.plot(img)


#step in time
##################################################
Nt = 12
times = np.linspace(0.0, 1.0/freq, Nt)*(solar_mass_per_s/mass)
phase = np.linspace(0.0, 1.0, Nt)

print 'angvel {}'.format(angvel)

for t, time_step in enumerate(times):
    print 'step: {:3d} / {:6.2f} / {:4.2f}'.format(t, time_step, phase[t])

    spot.star_time = time_step

    visz.star(img, spot)
    bounds = visz.spot_bounding_box()

    #integrate
    min_lims = [bounds[0], bounds[1]]
    max_lims = [bounds[2], bounds[3]]
            
    vals, errs = cubature( flux, 
                          2, 2, 
                          min_lims, max_lims,
                          relerr=1.0e-3,
                          maxEval=10000,
                          adaptive='p'
                          )
    
    print "vals", vals," ",errs

    pause(0.001)





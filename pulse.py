import sys
sys.path.append('/Users/natj/projects/arcmancer/lib/')
import pyarcmancer as pyac

from img import Imgplane
from lineprofile import *
from spot import Spot
from visualize import Visualize
import radiation

import numpy as np
import matplotlib as mpl
from pylab import *
from matplotlib import cm
import scipy.interpolate as interp
from cubature import cubature


from timeit import default_timer as timer

#from joblib import Parallel, delayed
#import multiprocessing

#num_cores = multiprocessing.cpu_count()
#print "num of cores {}", num_cores
#mpl.rcParams['image.cmap'] = 'inferno'



##################################################
# Set up figure & layout
fig = figure(figsize=(3,4), dpi=150)
mpl.rc('font', family='serif', size=7)
mpl.rc('xtick', labelsize='x-small')
mpl.rc('ytick', labelsize='x-small')
mpl.rcParams['image.cmap'] = 'inferno'



##################################################
# Star parameters
R = 15.0
M = 1.6
freq = 600.0
#freq = 1.0
incl = 60.0


##################################################
# Spot parameters
rho = 30.0
colat = 50.0


##################################################
#Physical conversion factors 
#G = 6.67384e-8
#c = 2.99792458e10
#Msun = 1.9885469e33

# Operate in units where G = c = 1.
# Use units of solar masses
solar_mass_per_km = 0.677220002407 # 1 / G*Msun/c^2
solar_mass_per_s  = 2.03e5


## conversion factors
kelvin_per_kev = 1.16045e7
km_per_kpc     = 3.08567758135e16
kg_per_Msun    = 1.988435e30
cm_tenkpc = 3.24077929e-23 #1cm/10kpc # 3.08567758135

# Variables in units of solar mass are derived here
# and presented with full name
mass        = M
radius      = R * solar_mass_per_km / mass
angvel      = freq * 2.0*np.pi / solar_mass_per_s * mass

imgscale    = (mass/solar_mass_per_km*1.0e5)**2 #cm^2/Msun



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


#Oblate Sch #WORKS
#metric = pyac.AGMMetric(radius, 1.0, angvel, pyac.AGMMetric.MetricType.agm_no_quadrupole)
#ns_surface = pyac.AGMSurface(radius, 1.0, angvel, pyac.AGMSurface.SurfaceType.agm_no_quadrupole)
#compactness = np.sqrt(1 - 2/radius) #isotropic radius compactness


#Full AGM metric & surface #WORKS
metric = pyac.AGMMetric(radius, 1.0, angvel, pyac.AGMMetric.MetricType.agm_standard)
ns_surface = pyac.AGMSurface(radius, 1.0, angvel, pyac.AGMSurface.SurfaceType.agm)
compactness = np.sqrt(1 - 2/radius) #isotropic radius compactness


surfaces = [ ns_surface ]



# Build and configure image plane by hand
img = Imgplane(conf, metric, surfaces)

img.verbose  = 1
img.incl     = np.deg2rad(incl) #set inclination
img.distance = 100.0*mass #set distance


#Locate star edges
img.find_boundaries(Nedge=10)

#Build internal coarse grid for the interpolation routines
img.generate_internal_grid(Nrad = 20, Nchi = 20 )
img.dissect_geos()


###################################################
# Radiation parameters
teff = 2.0 #Effective temperature
energies = np.array([2.0, 6.0, 12.0])  #Energy grid to compute the monochromatic fluxes
normalization = imgscale*cm_tenkpc**2 #normalization to physical cm^2 in the observers sky


###################################################
# flux function for integration
def flux(xy):
    x, y = xy

    time, phi, theta, cosa, reds = img.get_pixel(x, y) 
    #time, phi, theta, cosa, reds = img.get_exact_pixel(x, y) 

    EEd = reds #E_obs/E_surf

    coords = np.array([time, theta, phi])

    hit = spot.hit(coords)

    if hit:
        beam = radiation.isotropic_beaming(cosa)

        fluxNB = (EEd**3) * radiation.NB(teff) * beam
        fluxB  = (EEd**4) * radiation.EB(teff) * beam

        fluxE  = (EEd**3) * radiation.BE(teff, energies/EEd) * beam
        fluxNE = (EEd**2) * radiation.NBE(teff, energies/EEd) * beam

        fluxarray = np.hstack((fluxNE, fluxNB, fluxB, fluxE))*normalization

        return fluxarray
    else:
        return np.zeros(8)




#Construct output xy image plane from img object
##################################################
ion()
visz = Visualize()
visz.gs.update(hspace = 0.5)
visz.compactness = compactness
visz.plot(img)


# Set up pulse profile figure
visz.axs[5] = subplot( visz.gs[3,:] )
visz.axs[5].minorticks_on()
visz.axs[5].set_xlabel(r'Phase')
visz.axs[5].set_ylabel(r'Flux')
visz.axs[5].set_xlim(0,1)


#step in time
##################################################
Nt = 16
times = np.linspace(0.0, 1.0/freq, Nt)*(solar_mass_per_s/mass)
phase = np.linspace(0.0, 1.0, Nt)

fluxes = np.zeros((Nt, 9))


for t, time_step in enumerate(times):
    print 'step: {:3d} / {:6.2f} / {:4.2f}'.format(t, time_step, phase[t])

    spot.star_time = time_step

    visz.star(img, spot)

    visz.improve_on_spot_boundaries()
    bounds = visz.spot_bounding_box()

    #integrate
    min_lims = [bounds[0], bounds[1]]
    max_lims = [bounds[2], bounds[3]]

    start = timer()
    vals, errs = cubature( flux, 
                          2, 8, 
                          min_lims, max_lims,
                          relerr=1.0e-3,
                          maxEval=1000000,
                          adaptive='p'
                          )

    #Sometimes polynomial cubature does not find the spot so
    # we fall back to area splitted cubature. It is more robust
    # but slower & less accurate
    if max(vals) == 0.0:
        print "    fallback to area splitting cubature"
        vals, errs = cubature( flux, 
                              2, 8, 
                              min_lims, max_lims,
                              relerr=1.0e-3,
                              maxEval=1000000,
                              adaptive='h'
                              )
    end = timer()


    
    print '   flux {:6.2f} +- {:6.2f}%  | elapsed time: {}'.format(vals[3], 100.0*errs[3]/vals[3], end-start)


    #plot pulse profile on the fly
    fluxes[t,0]  = phase[t]
    fluxes[t,1::] = vals
    visz.axs[5].plot(phase[0:t], fluxes[0:t,4], "b.-")


    pause(0.001)


ioff()
#show()


#Finally save to file
fname = 'f{:03d}_bb_r{:02d}_m{:03.1f}_d{:02d}_i{:02d}_x{:02d}.csv'.format(
        np.int(freq),
        np.int(R),
        M,
        np.int(colat),
        np.int(incl),
        np.int(rho)
        )
print 'Saving to a file: '+fname


np.savetxt('out/'+fname,
        fluxes, 
        delimiter=',', 
        fmt = '%10.9e',
        header='Phase, F_N(2), F_N(6), F_N(12), F_Nb, F_N, F_E(2), F_E(6), F_E(12)'
        )


import sys
sys.path.append('/Users/natj/projects/arcmancer/lib/')
import pyarcmancer as pyac

from img import Imgplane
from visualize_polar import Visualize
from lineprofile import *
import units

import numpy as np
import matplotlib as mpl
from pylab import *
import os

from matplotlib import cm

import scipy.interpolate as interp

#from joblib import Parallel, delayed
#import multiprocessing


outdir = 'out/sweep/'


##################################################
# Set up figure & layout
fig = figure(figsize=(6,10)) 
mpl.rc('font', family='serif')
mpl.rc('xtick', labelsize='x-small')
mpl.rc('ytick', labelsize='x-small')
mpl.rcParams['image.cmap'] = 'inferno'


#num_cores = multiprocessing.cpu_count()
#print "num of cores {}", num_cores


#Setup pyarcmancer
##################################################
conf = pyac.Configuration()

conf.absolute_tolerance       = 1e-12
conf.relative_tolerance       = 1e-12
conf.henon_tolerance          = 1e-8
conf.sampling_interval        = 1e-3
conf.minimum_stepsize         = 1e-10
conf.maximum_steps            = 10000
conf.enforce_maximum_stepsize = False
conf.enforce_minimum_stepsize = True
conf.enforce_maximum_steps    = True
conf.store_only_endpoints     = True


#pyac.Log.set_console()
pyac.Log.set_file()



##################################################
# Star parameters
#R = 12.0
#M = 1.6
freq = 600.0
#incl = 15.0

for M in [1.5, 1.1, 1.8]:
    print "##################################################"
    print "M = ", M

    for R in [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0]:
        print "##################################################"
        print "  R = ", R
    
        #for incl in [90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 5, 1]:
        #for incl in [9, 8, 7, 6, 4, 3, 2, 0.5]:
        for incl in [90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 9, 8, 7, 6, 5, 4,3,2,1,0.5]:
            print "##################################################"
            print "    i = ",incl
        
            fname = 'neutronstar_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.png'.format( np.int(freq), np.int(R), M, np.int(incl))
            
            if os.path.isfile( outdir+fname ):
                continue

        
            # Variables in units of solar mass are derived here
            # and typically presented with full name
            mass        = M
            radius      = R * units.solar_mass_per_km / mass
            angvel      = freq * 2.0*np.pi / units.solar_mass_per_s * mass
            
            imgscale    = (mass/units.solar_mass_per_km*1.0e5)**2 #cm^2/Msun
            compactness = np.sqrt(1 - 2/radius) #isotropic radius compactness
            
            
            
            
            conf.absolute_tolerance       = 1e-12 * radius
            conf.minimum_stepsize         = 1e-10 * radius
            
            ##################################################
            #Define metric and surfaces of the spacetime 
            #S+D
            #metric = pyac.AGMMetric(radius, 1.0, angvel, pyac.AGMMetric.MetricType.agm_no_quadrupole)
            #ns_surface = pyac.AGMSurface(radius, 1.0, angvel, pyac.AGMSurface.SurfaceType.spherical)
            
            #Oblate Sch #WORKS
            #metric = pyac.AGMMetric(radius, 1.0, angvel, pyac.AGMMetric.MetricType.agm_no_quadrupole)
            #ns_surface = pyac.AGMSurface(radius, 1.0, angvel, pyac.AGMSurface.SurfaceType.agm_no_quadrupole)
            
            
            #Full AGM + oblate
            metric = pyac.AGMMetric(radius, 1.0, angvel, pyac.AGMMetric.MetricType.agm_standard)
            ns_surface = pyac.AGMSurface(radius, 1.0, angvel, pyac.AGMSurface.SurfaceType.agm)
            surfaces = [ ns_surface ]
            
            
            
            # Build and configure image plane by hand
            img = Imgplane(conf, metric, surfaces)
            
            img.verbose  = 1
            img.incl     = np.deg2rad(incl) #set inclination
            img.distance = 100000.0*mass #set distance
            
            
            #Locate star edges
            img.find_boundaries(Nedge=50, reltol=1.0e-4, max_iterations=30)
            
            #Build internal coarse grid for the interpolation routines
            img.generate_internal_grid(Nrad = 80, Nchi = 50 )
            img.dissect_geos()
            
            
            
            #Construct output xy image plane from img object
            ##################################################
            ion()
            visz = Visualize()
            visz.gs.update(hspace = 0.5)
            visz.compactness = compactness
            visz.plot(img)

            
            #prepare line profile axis object
            visz.axs[6] = subplot( visz.gs[3,:] )
            visz.axs[6].minorticks_on()
            visz.axs[6].set_xlabel(r'Energy')
            visz.axs[6].set_ylabel(r'Flux')
            
            
            
            #Construct image
            #visz.star(img, spot)
            #visz.polar(img, spot)
            
            visz.dissect(img)
            visz.star_plot(0.0)
            
            visz.polar_dissect(img)
            visz.polar_plot(0.0)
            
            
            
            ##################################################
            # Compute line profile
            es, yy2 = lineprofile(visz.redshift**4, visz.redshift)
            dE = np.max( np.abs(es[0] - compactness), np.abs(compactness - es[-1]))
            
            
            
            ##################################################
            #Save redshift into a file
            fname = 'reds_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.csv'.format(
                    np.int(freq),
                    np.int(R),
                    M,
                    np.int(incl),
                    )
            print 'Saving to a file: '+fname
            
            np.savetxt(outdir+fname,
                    visz.redshift.flatten(),
                    delimiter=',', 
                    fmt = '%10.9e'
                    )
            
            #Save thetas into a file
            fname = 'thetas_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.csv'.format(
                    np.int(freq),
                    np.int(R),
                    M,
                    np.int(incl),
                    )
            print 'Saving to a file: '+fname
            
            np.savetxt(outdir+fname,
                    visz.thetas.flatten(),
                    delimiter=',', 
                    fmt = '%10.9e'
                    )
            
            #Save phi into a file
            fname = 'phis_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.csv'.format(
                    np.int(freq),
                    np.int(R),
                    M,
                    np.int(incl),
                    )
            print 'Saving to a file: '+fname
            
            np.savetxt(outdir+fname,
                    visz.phis.flatten(),
                    delimiter=',', 
                    fmt = '%10.9e'
                    )
            
            
            #redshift limits
            vmin = compactness - dE
            vmax = compactness + dE
            
            
            # Line profile
            ##################################################
            #ax = subplot(gs[2,2])
            #ax.set_xlim(0.8, 1.2)
            visz.axs[6].plot(es, yy2, "b-")
            
            
            pause(1.0)
            fname = 'neutronstar_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.png'.format(
                    np.int(freq),
                    np.int(R),
                    M,
                    np.int(incl),
                    )
            
            savefig(outdir+fname)
            
            
            #save lineprofile
            ##################################################
            
            
            #Finally save to file
            fname = 'lineprofile_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.csv'.format(
                    np.int(freq),
                    np.int(R),
                    M,
                    np.int(incl),
                    )
            print 'Saving to a file: '+fname
            
            
            np.savetxt(outdir+fname,
                    np.vstack((es, yy2)).T, 
                    delimiter=',', 
                    fmt = '%10.9e',
                    header='Energy, pdf'
                    )
        
        

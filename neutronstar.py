import sys
sys.path.append('/Users/natj/projects/arcmancer/lib/')
import pyarcmancer as pyac

from img import Imgplane
from lineprofile import *

import numpy as np
import matplotlib as mpl
from pylab import *

from matplotlib import cm

import scipy.interpolate as interp

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
M = 1.6
freq = 400.0
incl = 15.0



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
img.find_boundaries(Nedge=100, reltol=1.0e-4, max_iterations=30)

#Build internal coarse grid for the interpolation routines
img.generate_internal_grid(Nrad = 100, Nchi = 100 )

img.dissect_geos()



#Construct output xy image plane from img object
##################################################
#x_span = 1.5*radius
#y_span = 1.5*radius
x_span = 11.0
y_span = 11.0


x_bins = 500
y_bins = 500

pixel_dx = 2*x_span / x_bins
pixel_dy = 2*y_span / y_bins
pixel_area = pixel_dx * pixel_dy


xs = np.linspace(-x_span, x_span, x_bins)
ys = np.linspace(-y_span, y_span, y_bins)


#build interpolated array
redshift      = np.zeros((x_bins, y_bins))
obs_hit_angle = np.zeros((x_bins, y_bins))
times         = np.zeros((x_bins, y_bins))
thetas        = np.zeros((x_bins, y_bins))
phis          = np.zeros((x_bins, y_bins))



# Finally, loop over pixels to build image
for i, xi in enumerate(xs):
    if i % 10 == 0:
        print "{} % done".format(float(i)/len(xs) * 100)

    for j, yi in enumerate(ys):
        time, phi, theta, cosa, reds = img.get_pixel(xi, yi)


        #time, phi, theta, cosa, reds = img.get_exact_pixel(xi, yi)

        redshift[i,j]      = reds
        obs_hit_angle[i,j] = cosa
        times[i,j]         = time
        thetas[i,j]        = theta
        phis[i,j]          = phi



##################################################
# Compute line profile
es, yy2 = lineprofile(redshift**4, redshift)
#es = es/compactness
dE = np.max( np.abs(es[0] - compactness), np.abs(compactness - es[-1]))


##################################################
# plot values on image plane
def trans(mat):
    return np.flipud(mat.T)
#return mat
def detrans(mat):
    return np.flipud(mat).T

def clean_image(mat):

    #mask all 0.0 elements and transpose
    mat_masked = np.ma.masked_where(mat == 0, mat) 

    return trans(mat_masked)

#build up a chess board layering from phi and theta coordinates
def chess_layer(phis, thetas):

    none  = 0
    white = 1
    black = 2
    
    Nx, Ny = np.shape(phis)
    mat    = np.zeros((Nx,Ny))

    for i in range(Nx):
        for j in range(Ny):
            phi = phis[i,j]
            theta = thetas[i,j]
            
            if (phi == 0.0) and (theta == 0):
                mat[i,j] = none
                continue

            xd = np.int(60.0*phi/(2*pi))
            yd = np.int(60.0*theta/(2*pi))
            
            if (xd & 1) and (yd & 1):
                mat[i,j] = black
            elif (xd & 0) and (yd & 0):
                mat[i,j] = black
            else:
                mat[i,j] = white

    return mat

#Compute chess pattern
chess = chess_layer(phis, thetas)


##################################################
#Save redshift into a file
fname = 'reds_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.csv'.format(
        np.int(freq),
        np.int(R),
        M,
        np.int(incl),
        )
print 'Saving to a file: '+fname

np.savetxt('out/'+fname,
        redshift.flatten(),
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

np.savetxt('out/'+fname,
        thetas.flatten(),
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

np.savetxt('out/'+fname,
        phis.flatten(),
        delimiter=',', 
        fmt = '%10.9e'
        )




#Clean images for imshow
obs_hit_angle = clean_image(obs_hit_angle)
redshift      = clean_image(redshift)
times         = clean_image(times)
thetas        = clean_image(thetas)
phis          = clean_image(phis)
chess         = clean_image(chess)


# other settings for imshow
extent=( xs[0], xs[-1], ys[0], xs[-1] )
interpolation = 'nearest'

#redshift limits
vmin = compactness - dE
vmax = compactness + dE

# General image
ax = subplot(gs[0:2,0:2])
ax.axis('off')
ax.imshow(chess, interpolation=interpolation, extent=extent, cmap=cm.get_cmap('Greys'), vmin=0.8, vmax=2.0, alpha=0.6)
ax.imshow(redshift, interpolation=interpolation, origin='lower', extent=extent,
        cmap=cm.get_cmap('coolwarm_r'), vmin=vmin, vmax=vmax, alpha=0.95)

levels = np.linspace(vmin, vmax, 20)
ax.contour(redshift, levels, hold='on', colors='w',
        origin='lower', extent=extent, vmin=0.8*compactness, vmax=1.2*compactness)



# Observer hit angle (cos\alpha)
ax = subplot(gs[2,0])
ax.minorticks_on()
cax = ax.imshow(obs_hit_angle, interpolation=interpolation, extent=extent)
colorbar(cax)
ax.set_title(r'emitter angle $\alpha$')


# Redshift with contours
ax = subplot(gs[2,1])
ax.minorticks_on()
#cax = ax.imshow(redshift, interpolation=interpolation, origin='lower', extent=extent,
#        cmap=cm.get_cmap('coolwarm_r'), vmin=vmin, vmax=vmax)
#
#levels = np.linspace(vmin, vmax, 20)
#ax.contour(redshift, levels, hold='on', colors='w',
#        origin='lower', extent=extent, vmin=vmin, vmax=vmax)

# Redshift^4 (monoenergetic flux)
cax = ax.imshow(redshift**4, interpolation=interpolation, origin='lower', extent=extent,
        cmap=cm.get_cmap('coolwarm_r'))
ax.contour(redshift**4, 20, hold='on', colors='w',
        origin='lower', extent=extent)

colorbar(cax)
ax.set_title('redshift')

# Phi angle
ax = subplot(gs[0,2])
ax.minorticks_on()
cax = ax.imshow(phis, interpolation=interpolation, extent=extent)
colorbar(cax)
ax.set_title(r'$\phi$')


# Theta angle
ax = subplot(gs[1,2])
ax.minorticks_on()
cax = ax.imshow(thetas, interpolation=interpolation, extent=extent)
colorbar(cax)
ax.set_title(r'$\theta$')


# Line profile
ax = subplot(gs[2,2])
#ax.set_xlim(0.8, 1.2)
ax.plot(es, yy2, "b-")
ax.set_title(r'line profile')


savefig('neutronstar.png')


#save lineprofile


#Finally save to file
fname = 'lineprofile_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.csv'.format(
        np.int(freq),
        np.int(R),
        M,
        np.int(incl),
        )
print 'Saving to a file: '+fname


np.savetxt('out/'+fname,
        np.vstack((es, yy2)).T, 
        delimiter=',', 
        fmt = '%10.9e',
        header='Energy, pdf'
        )



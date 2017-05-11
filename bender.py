import sys
sys.path.append('/Users/natj/projects/arcmancer/lib/')
import pyarcmancer as pyac
import numpy as np

import matplotlib as mpl
#import matplotlib.pyplot as plt
from pylab import *

from matplotlib import cm

#from scipy import interpolate
import scipy.interpolate as interp

from joblib import Parallel, delayed
import multiprocessing


mpl.rc('font', family='serif')
mpl.rc('xtick', labelsize='small')
mpl.rc('ytick', labelsize='small')
gs = GridSpec(3, 3)
gs.update(hspace = 0.3)
#gs.update(wspace = 0.3)

num_cores = multiprocessing.cpu_count()
print "num of cores {}", num_cores

mpl.rcParams['image.cmap'] = 'inferno'





##################################################
#Setup star
R_in            = 12.0
m_in            = 1.4
freq_in         = 600.0
colat_in        = 10.0

spot_temp_in    = 2.0
spot_ang_in     = 10.0
spot_colat_in   = 90.0




##################################################
#Planck units
## physical constants
## since c = G = 1, we get c^2/G = 0.677 Msun / km
L_to_M = 0.677

## c in kilometers per second
c = 299792.458

## G in km^3 kg^-1 second^-2
G = 1.327e11

## Stefan's constant in keV / (s cm^2 K^4)
sigma_sb = 35392.0

#pi from numpy
pi = np.pi

## conversion factors
kelvin_per_kev = 1.16045e7
km_per_kpc     = 3.086e16
kg_per_Msun    = 1.988435e30


# Operate in units where G = c = 1.
# Use units of solar masses
solar_mass_per_km = 0.6772
solar_mass_per_s  = 2.03e5


#G=c=1 units (i.e., solar masses)
mass = m_in
R_eq = R_in * solar_mass_per_km / mass
angvel = freq_in * 2.0*np.pi / solar_mass_per_s * mass

compactness = np.sqrt(1 - 2/R_eq)
print "compactness=",compactness

##################################################
#convert spot quantities
spot_temp  = spot_temp_in * kelvin_per_kev
spot_colat = spot_colat_in * pi/180.0
spot_ang   = spot_ang_in * pi/180.0

#define luminosity distance of NS = 10 kpc
#luminosity_distance = 10.0 * km_per_kpc * 1e3 / planck_length_in_m


#Setup pyarcmancer
##################################################
conf = pyac.Configuration()

conf.absolute_tolerance = 1e-12 * R_eq
conf.relative_tolerance = 1e-12
conf.henon_tolerance    = 1e-8
conf.sampling_interval = 1e-3
conf.minimum_stepsize  = 1e-10 * R_eq
conf.maximum_steps     = 10000
conf.enforce_maximum_stepsize = False
conf.enforce_minimum_stepsize = True
conf.enforce_maximum_steps    = True
conf.store_only_endpoints     = True



#Define metric of the spacetime
#metric = pyac.SchwarzschildMetric(mass/R_eq)
metric = pyac.AGMMetric(R_eq, 1.0, angvel, pyac.AGMMetric.MetricType.agm_standard)


ns_surface = pyac.AGMSurface(R_eq, 1.0, angvel, pyac.AGMSurface.SurfaceType.spherical)
#ns_surface = pyac.AGMSurface(R_eq, 1.0, angvel, pyac.AGMSurface.SurfaceType.agm)
surfaces = [ ns_surface ]



#pyac.Log.set_console()
pyac.Log.set_file()


##################################################
# Construct image plane
inclination = np.deg2rad(colat_in)
distance = mass * 100

x_span = 1.5*R_eq
y_span = 1.5*R_eq

x_bins = 200
y_bins = 200

pixel_dx = 2*x_span / x_bins
pixel_dy = 2*y_span / y_bins
pixel_area = pixel_dx * pixel_dy


initial_xs = np.linspace(-x_span, x_span, x_bins)
initial_ys = np.linspace(-y_span, y_span, y_bins)


# construct local spherical axes in cartesian coordinates
def local_spherical_axes(pos_cart):
    #print "Computing cartesian components of local spherical axes at {}".format(pos_cart)
    u_r = pos_cart / np.linalg.norm(pos_cart)
    u_phi = np.cross(np.array([0,0,1]), u_r)
    u_phi /= np.linalg.norm(u_phi)
    u_theta = -np.cross(u_r, u_phi)
    u_theta /= np.linalg.norm(u_theta)
    #print "result {} {} {}".format(u_r, u_theta, u_phi)
    return [u_r, u_theta, u_phi]


# in the limit m -> 0, BL coordinates go to oblate spheroidal minkowski. These
# go to minkowski for r -> inf or a ->0
def boyer_lindquist_position(x_cart):
    #print "Transforming cartesian position {} to Boyer-Lindquist".format(x_cart)
    x, y, z = x_cart
    r = np.linalg.norm(x_cart)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    return np.array([0, r, theta, phi])


# in the limit m = 0
def cartesian_position(a, x_bl):
    r, theta, phi = x_bl[1:]

    x = np.hypot(r, a) * np.sin(theta) * np.cos(phi)
    y = np.hypot(r, a) * np.sin(theta) * np.sin(phi)
    z = r*np.cos(theta)
    return np.array([x_bl[0], x, y, z])

# Initialize photon with some (x,y) coordinates in the _image plane_
# and make it point towards the neutron star
def xy2geo(metric, distance, inclination, x, y):

    #print "Initializing geodesic for {},{}".format(x,y)

    # get coordinates for position of image plane point
    normal = np.array([np.sin(inclination), 0.0, np.cos(inclination)])

    x_cart = np.array([-y*np.cos(inclination), x, y*np.sin(inclination)])
    #x_cart = np.array([ x, -y*np.cos(inclination), y*np.sin(inclination) ] )



    x_cart += distance * normal
    x_sph = boyer_lindquist_position(x_cart)


    # get velocity by projecting to local B-L axes
    # get axes in _cartesian_ coordinates
    u_r, u_theta, u_phi = local_spherical_axes(x_cart)
    vel_cart = normal
    vel_sph = np.array([0, 
                        np.dot(u_r     , vel_cart) ,
                        np.dot(u_theta , vel_cart) / x_sph[1],
                        np.dot(u_phi   , vel_cart) / (x_sph[1] * np.sin(x_sph[2]))])


    # define vertical and horizontal
    vert = pyac.normalize(metric, x_sph, np.array([0, 0, -1.0, 0]))

    vert_vel = pyac.project_along(metric, x_sph, vel_sph, vert)
    vert -= vert_vel
    vert = pyac.normalize(metric, x_sph, vert)

    horz = pyac.spatial_cross_product(
        metric, pyac.static_observer(metric, x_sph), vert, vel_sph)
    horz = pyac.normalize(metric, x_sph, horz)
    horz = pyac.normalize(metric, x_sph, np.array([0, 0, 0, 1.0]))

    if 0:
        # test
        print "vert", vert, "horz", horz
        print "vert.horz", metric.dot(x_sph, vert, horz)
        print "vert.u", metric.dot(x_sph, vert, vel_sph)
        print "horz.u", metric.dot(x_sph, horz, vel_sph)


    geo = pyac.Geodesic(metric, x_sph, vel_sph, vert, horz, pyac.VectorType.null)

    return geo


#polar coordinates to photon in image plane
def pol2geo(metric, distance, inclination, rad, chi):
    x = rad * np.sin(chi)
    y = rad * np.cos(chi)

    return xy2geo(metric, distance, inclination, x, y)


# generate an image plane of geodesics
def generate_image_plane(metric, distance, inclination, x_span, y_span, x_bins, y_bins):
    xs = initial_xs
    ys = initial_ys

    # avoid line of poles?
    #xs += 0.1*pixel_dx

    plane = []
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            #print "ix", ix, "x", x
            #print "iy", iy, "y", y

            #Make photon from (x,y) coords
            geo = xy2geo(metric, distance, inclination, x, y)
            plane.append((ix, iy, x, y, geo))

    return plane


def compute_element(el, distance, metric, conf, surfaces):
    el[4].compute(-(distance + 5.0*R_eq), metric, conf, surfaces)


# black-body specific intensity in natural units (h = 1)
def bb_intensity(nu, T):
    return 2 * nu**3 * ( np.exp( nu/T ) - 1 )**(-1)



class ImagePlane:

    def __init__(self, metric, dist, incl):

        self.metric = metric
        self.dist   = dist
        self.incl   = incl


class Pixel:

    ix = 0
    iy = 0

    def __init__(self, x, y):

        self.x  = x
        self.y  = y
        
        #self.geo = geo


##################################################

#Find star radius boundaries
def find_boundaries(metric, distance, inclination, surfaces):
    print "Finding edge boundaries for the star..."

    Nedge = 10
    chis = np.linspace(0.0, 1.0, Nedge)*2.0*pi + 0.001
    rlims = np.zeros(Nedge)

    rmin = 0.0
    rmax = 12.0


    for i, chii in enumerate(chis):
        geos = []

        rmini = rmin
        rmaxi = rmax
        rmid = 0.0

        Nbi = 20
        N = 0
        #for N in range(Nbi):

        relerr = 1.0
        reltol = 1e-3
        rmid_old = 100.0

        while (N < Nbi) and (relerr > reltol):
            rmid = (rmini + rmaxi)/2.0

            geos.append((0,0,0,0, pol2geo(metric, distance, inclination, rmid, chii) ))
            compute_element(geos[N], distance, metric, conf, surfaces)


            hit = geos[N][4].front_termination().hit_surface

            if hit:
                rmini = rmid
            else:
                rmaxi = rmid
            
            relerr = np.abs(rmid - rmid_old)/rmid
            rmid_old = rmid

            N += 1
            print "Iterating edge at {} after {} tries for angle={} ({})".format(rmid, N, chii, relerr)

        rlims[i] = rmid
    return chis, rlims
                

##################################################
# Creates internal polar grid for the interpolation
def internal_polar_grid(Nrad, Nchi):
    
    dchi_edge = 0.001
    chimin = 0.0 - dchi_edge
    chimax = 2.0*pi + dchi_edge
    
    chi_diffs = 0.8 + np.sin( np.linspace(0.0, 2.0*pi, Nchi-3) )**2
    chi_diffs = np.insert(chi_diffs, 0, 0.0)
    chi_grid = chimin + (chimax - chimin) * np.cumsum(chi_diffs)/np.sum(chi_diffs)
    chi_grid = np.insert(chi_grid, 0, chi_grid[0] - dchi_edge)
    chi_grid = np.append(chi_grid, chi_grid[-1] + dchi_edge)

    #chi_grid = np.linspace(0.0, 2.0*pi, Nchi)
    
    
    rad_diffs = 1.0 / np.exp( np.linspace(1.2, 2.0, Nrad-1)**2)
    rad_grid = rmax * np.cumsum(rad_diffs)/np.sum(rad_diffs)
    rad_grid = np.insert(rad_grid, 0, 0.001)
    
    
    grid = np.empty((Nrad,Nchi), dtype=np.object)
    
    for i, chi in enumerate(chi_grid):
        print "{} % done".format(float(i)/len(chi_grid) * 100)
        for j, rad in enumerate(rad_grid):
            print "  tracing geodesic at chi={} and r={}".format(chi, rad) 
            #Trace geodesic from image plane to star
            grid[j,i] = pol2geo(metric, distance, inclination, rad, chi)
            grid[j,i].compute(-(distance + 5.0*R_eq), metric, conf, surfaces)
    
    return rad_grid, chi_grid, grid


#Reduce geodesic into auxiliary quantities
def dissect_geos(grid, rad_grid, chi_grid):

    print "Dissecting geodesic paths to observed quantities..."

    Nrad, Nchi = np.shape(grid)

    Reds   = np.zeros((Nrad, Nchi))
    Cosas  = np.zeros((Nrad, Nchi))
    Times  = np.zeros((Nrad, Nchi))
    Thetas = np.zeros((Nrad, Nchi))
    Phis   = np.zeros((Nrad, Nchi))

    
    for i, chi in enumerate(chi_grid):
        #print "{} % done".format(float(i)/len(chi_grid) * 100)
        for j, rad in enumerate(rad_grid):

            geo = grid[j,i]

            hit_pt = geo.get_points()[0]
            obs_pt = geo.get_points()[-1]

            #coordinates 
            #surface_point.point.x
            t  = Times[j,i]  = hit_pt.point.x[0]
            th = Thetas[j,i] = hit_pt.point.x[2]
            p  = Phis[j,i]   = hit_pt.point.x[3]
            
            hit = geo.front_termination().hit_surface

            #Redshift
            g = Reds[j,i] = \
                metric.dot(obs_pt.point, pyac.static_observer(metric, obs_pt.x())) / \
                    metric.dot(hit_pt.point, ns_surface.observer(metric, hit_pt.x()))

            #hit angle
            cosa = Cosas[j,i] = \
                geo.front_termination().observer_hit_angle


    return Reds, Cosas, Times, Thetas, Phis


def mod2pi(x):
    #return np.fmod(x, 2*np.pi)

    while (x > 2.0*np.pi):
        x -= 2.0*np.pi
    while (x < 0.0):
        x += 2.0*np.pi
    return x



def calc_chi(x,y):
    return mod2pi(np.pi/2.0 - np.arctan2(y,x) )
    #return mod2pi( np.arctan2(y,x) )


#Interpolation function
# We take x and y grid values in polar coordinates
# New points (that we are interpolating into) are
# in Cartesian (x,y) coordinates.
# New points are then first transformed to polar coordinates for querying.
def grid_interpolate(rads, chis, z, 
                     x2, y2,
                     edge_inter):

    ##rads = np.hypot(x, y)
    ##chis = np.arctan2(y, x)

    ##rads2 = np.hypot(x2, y2)
    ##chis2 = np.arctan2(y2, x2)

    #x_sparse, y_sparse = np.meshgrid(rads, chis)
    ##y_sparse, x_sparse = np.meshgrid(rads, chis)

    ##x_dense, y_dense = np.meshgrid(rads2, chis2)

    ##new way of transforming only after the unpacking of mesh
    #x_dense, y_dense = np.meshgrid(x2, y2)


    #z2 = interp.griddata(np.array([x_sparse.ravel(),y_sparse.ravel()]).T,
    #                              z.ravel(),
    #                              (x_dense,y_dense), method='linear')
    #print np.shape(z)
    #print np.shape(rads), min(rads), max(rads)
    #print np.shape(chis), min(chis), max(chis)

    #Using spline; XXX is this a good idea, I dunno!
    ir = interp.RectBivariateSpline(rads, chis, z, 
                                    #bbox=[0.01, max(rads), 0.0, 2*np.pi], 
                                    kx=1, ky=1, s=0)

    #build interpolated array
    z2 = np.zeros((len(x2), len(y2)))
    for i, xi in enumerate(x2):
        for j, yi in enumerate(y2):
            radi = np.hypot(xi, yi)
            #chii = np.arctan2(yi, xi)
            chii = calc_chi(xi, yi)

            #Check radius
            redge = edge_inter(chii)
            if radi <= redge:
                z2[i,j] = ir.ev(radi, chii)
                #z2[i,j] = chii


    return z2

##################################################
##################################################

#Locate edges
chis, rlims = find_boundaries(metric, distance, inclination, surfaces)
#rlims = 8.04
rmax = np.max(rlims)*1.001
print "Maximum edge {}".format(rmax)

#Build edge location interpolator
edge_inter_raw = interp.InterpolatedUnivariateSpline(chis, rlims)

#Build internal coarse grid for the interpolation routines
rad_grid, chi_grid, coarse_polar_grid = internal_polar_grid(30, 30)

#Read observed values from geodesics
reds_int, cosas_int, times_int, thetas_int, phis_int  = dissect_geos(coarse_polar_grid, rad_grid, chi_grid)

#Interpolate everything into thick grid
redshift = grid_interpolate(rad_grid, chi_grid, reds_int,
                            initial_xs, initial_ys, edge_inter_raw)
obs_hit_angle = grid_interpolate(rad_grid, chi_grid, cosas_int,
                            initial_xs, initial_ys, edge_inter_raw)
times    = grid_interpolate(rad_grid, chi_grid, times_int,
                            initial_xs, initial_ys, edge_inter_raw)
thetas   = grid_interpolate(rad_grid, chi_grid, thetas_int,
                            initial_xs, initial_ys, edge_inter_raw)
phis     = grid_interpolate(rad_grid, chi_grid, phis_int,
                            initial_xs, initial_ys, edge_inter_raw)


# Computes lineprofile
def lineprofile(fluxc, redsc):
    print "Computing line profile..."

    xarr = np.array([])
    yarr = np.array([])

    Ny_dense, Nx_dense = np.shape(fluxc)

    energ = 1.0
    for jj in range(Ny_dense):
        for ii in range(Nx_dense):
            fy  = fluxc[jj, ii]
            xii = redsc[jj, ii]

            if xii > 0.0:
                xarr = np.append(xarr, xii*energ)
                yarr = np.append(yarr, fy)


    xind = np.argsort(xarr)
    xarrs = xarr[xind]
    yarrs = yarr[xind]

    NN = len(xarrs)

    emin = np.min(xarrs)*0.99
    emax = np.max(xarrs)*1.01

    Nr = 100
    es = np.linspace(emin, emax, Nr)
    yy2 = np.zeros((Nr))

    xst = 0
    for ii in range(1,Nr):
        for jj in range(xst, NN):
            if es[ii-1] <= xarrs[jj] < es[ii]:
                yy2[ii] += yarrs[jj]
            elif xarrs[jj] >= es[ii]:
                xst = jj
                break

    #normalize
    des = np.diff(es)[1]
    yy2 = yy2 / np.sum(yy2*des)

        
    return es, yy2         
 
es, yy2 = lineprofile(redshift**3, redshift)
es = es/compactness


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
    #return mat_masked

#build up a chess board layering from phi and theta coordinates
def chess_layer(phis, thetas):

    none = 0
    white = 1
    black = 2
    
    Nx, Ny = np.shape(phis)
    mat = np.zeros((Nx,Ny))

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


# transform everything to mesh with imshow

chess         = chess_layer(phis, thetas)

obs_hit_angle = clean_image(obs_hit_angle)
redshift      = clean_image(redshift)
times         = clean_image(times)
thetas        = clean_image(thetas)
phis          = clean_image(phis)
chess         = clean_image(chess)


# other settings for imshow
extent=( initial_xs[0], initial_xs[-1], initial_ys[0], initial_ys[-1])
interpolation = 'nearest'


ax = subplot(gs[0:2,0:2])
ax.axis('off')
ax.imshow(chess, interpolation=interpolation, extent=extent, cmap=cm.get_cmap('Greys'), vmin=0.8, vmax=2.0, alpha=0.6)
ax.imshow(redshift, interpolation=interpolation, origin='lower', extent=extent,
        cmap=cm.get_cmap('coolwarm_r'), vmin=0.8*compactness, vmax=1.2*compactness, alpha=0.95)
ax.contour(redshift, 10, hold='on', colors='w',
        origin='lower', extent=extent, vmin=0.8*compactness, vmax=1.2*compactness)



ax = subplot(gs[2,0])
ax.minorticks_on()
cax = ax.imshow(obs_hit_angle, interpolation=interpolation, extent=extent)
colorbar(cax)
ax.set_title(r'emitter angle $\alpha$')


ax = subplot(gs[2,1])
ax.minorticks_on()
cax = ax.imshow(redshift, interpolation=interpolation, origin='lower', extent=extent,
        cmap=cm.get_cmap('coolwarm_r'), vmin=0.8*compactness, vmax=1.2*compactness)
ax.contour(redshift, 20, hold='on', colors='w',
        origin='lower', extent=extent, vmin=0.8*compactness, vmax=1.2*compactness)
colorbar(cax)
ax.set_title('redshift')


ax = subplot(gs[0,2])
ax.minorticks_on()
cax = ax.imshow(phis, interpolation=interpolation, extent=extent)
colorbar(cax)
ax.set_title(r'$\phi$')


ax = subplot(gs[1,2])
ax.minorticks_on()
cax = ax.imshow(thetas, interpolation=interpolation, extent=extent)
colorbar(cax)
ax.set_title(r'$\theta$')


ax = subplot(gs[2,2])
ax.plot(es, yy2, "b-")
ax.set_title(r'line profile')


#plt.tight_layout()
#plt.show()
savefig('arcmancer_debug.png')

import sys
sys.path.append('/Users/natj/projects/arcmancer/lib/')
import pyarcmancer as pyac
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy import interpolate

from joblib import Parallel, delayed
import multiprocessing


num_cores = multiprocessing.cpu_count()
print "num of cores {}", num_cores

mpl.rcParams['image.cmap'] = 'inferno'





##################################################
#Setup star
R_in            = 12.0
m_in            = 1.4
freq_in         = 0.0
colat_in        = 90.0

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

## planck units
#planck_length_in_m = 1.616e-35
#planck_mass_in_kg  = 2.176e-8
#planck_time_in_s   = 5.3912e-44
#planck_temp_in_K   = 1.417e32
#planck_charge_in_C = 1.876e-18

# Operate in units where G = c = 1.
# Use units of solar masses
solar_mass_per_km = 0.6772
solar_mass_per_s  = 2.03e5


#G=c=1 units (i.e., solar masses)
mass = m_in
R_eq = R_in * solar_mass_per_km / mass
angvel = freq_in * 2.0*np.pi / solar_mass_per_s * mass


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
conf.maximum_steps     = 100000
conf.enforce_maximum_stepsize = False
conf.enforce_minimum_stepsize = True
conf.enforce_maximum_steps    = True
conf.store_only_endpoints     = False



#Define metric of the spacetime
#metric = pyac.SchwarzschildMetric(mass/R_eq)
metric = pyac.AGMMetric(R_eq, 1.0, angvel, pyac.AGMMetric.MetricType.agm_standard)


ns_surface = pyac.AGMSurface(R_eq, 1.0, angvel, pyac.AGMSurface.SurfaceType.spherical)
surfaces = [ ns_surface ]



#pyac.Log.set_console()
pyac.Log.set_file()


##################################################
# Construct image plane
inclination = np.deg2rad(colat_in)
distance = mass * 100

x_span = 1.5*R_eq
y_span = 1.5*R_eq

x_bins = 20
y_bins = 20

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
    Nedge = 5
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
        for N in range(Nbi):
            rmid = (rmini + rmaxi)/2.0

            geos.append((0,0,0,0, pol2geo(metric, distance, inclination, rmid, chii) ))
            compute_element(geos[N], distance, metric, conf, surfaces)


            hit = geos[N][4].front_termination().hit_surface

            if hit:
                rmini = rmid
            else:
                rmaxi = rmid
            #print "Iterating edge at {} after {} tries for {}".format(rmid, N, chii)

        rlims[i] = rmid
    return chis, rlims
                

chis, rlims = find_boundaries(metric, distance, inclination, surfaces)
rmax = np.max(rlims)*1.001
print "Maximum edge {}".format(rmax)


edge_inter_raw = interpolate.InterpolatedUnivariateSpline(chis, rlims)
#def edge_inter



##################################################
# Now create internal polar grid
Nchi = 20
Nrad = 20

dchi_edge = 0.001
chimin = 0.0 - dchi_edge
chimax = 2.0*pi + dchi_edge

chi_diffs = 0.8 + np.sin( np.linspace(0.0, 2.0*pi, Nchi-3) )**2
chi_diffs = np.insert(chi_diffs, 0, 0.0)
chi_grid = chimin + (chimax - chimin) * np.cumsum(chi_diffs)/np.sum(chi_diffs)
chi_grid = np.insert(chi_grid, 0, chi_grid[0] - dchi_edge)
chi_grid = np.append(chi_grid, chi_grid[-1] + dchi_edge)

print chi_grid


rad_diffs = 1.0 / np.exp( np.linspace(1.2, 2.0, Nrad-1)**2)
rad_grid = rmax * np.cumsum(rad_diffs)/np.sum(rad_diffs)
rad_grid = np.insert(rad_grid, 0, 0.001)

print rad_grid




##################################################
imgplane = generate_image_plane(metric, distance, inclination, x_span,
                                y_span, x_bins, y_bins)


#ns_surface = pyac.AGMSurface(R_eq, 1.0, angvel, pyac.AGMSurface.SurfaceType.spherical)
#surfaces = [ ns_surface ]


for i, el in enumerate(imgplane):
    if 10*i % len(imgplane) == 0:
        print "{} % done".format(float(i)/len(imgplane) * 100)
    compute_element(el, distance, metric, conf, surfaces)


# plot values on image plane
def trans(mat):
    return np.flipud(mat.T)
#return mat
def detrans(mat):
    return np.flipud(mat).T

xs             = np.zeros((x_bins, y_bins))
ys             = np.zeros((x_bins, y_bins))
iis            = np.zeros((x_bins, y_bins))
jjs            = np.zeros((x_bins, y_bins))
radii          = np.zeros((x_bins, y_bins))
obs_hit_angle  = np.zeros((x_bins, y_bins))
rest_hit_angle = np.zeros((x_bins, y_bins))
redshift       = np.zeros((x_bins, y_bins))
flux           = np.zeros((x_bins, y_bins))
bolflux        = np.zeros((x_bins, y_bins))
polfrac        = np.zeros((x_bins, y_bins))
polangle       = np.zeros((x_bins, y_bins))


for el in imgplane:
    i, j = el[0:2]
    hit_pt = el[4].get_points()[0]
    radii[i,j] = hit_pt.x()[1]
    iis[i,j] = i
    jjs[i,j] = j
    # _screen_ x and y
    xs[i,j] = el[2]
    ys[i,j] = el[3]
    #print "i {} j {} x {} y {}".format(i, j, xs[i,j], ys[i,j])


for el in imgplane:
    i, j = el[0:2]

    if not el[4].front_termination().hit_surface:
        #print "ray at {},{}", i, j, "didn't hit anything!"
        #print el[4].get_points()
        continue

    hit_pt = el[4].get_points()[0]
    obs_pt = el[4].get_points()[-1]


    obs_hit_angle[i,j] = el[4].front_termination().observer_hit_angle
    rest_hit_angle[i,j] = el[4].front_termination().rest_hit_angle


    g = redshift[i,j] = \
        metric.dot(obs_pt.point, pyac.static_observer(metric, obs_pt.x())) / \
        metric.dot(hit_pt.point, ns_surface.observer(metric, hit_pt.x()))


    #bolflux[i,j] = g**4 * ntdisk['flux'][i,j] # * pixel_area
    bolflux[i,j] = g**4 * 1.0


# transform everything to mesh with imshow
xs             = trans(xs            )
ys             = trans(ys            )
iis            = trans(iis           )
jjs            = trans(jjs           )
radii          = trans(radii         )
obs_hit_angle  = trans(obs_hit_angle )
rest_hit_angle = trans(rest_hit_angle)
redshift       = trans(redshift      )
flux           = trans(flux          )
bolflux        = trans(bolflux       )
polfrac        = trans(polfrac       )
polangle       = trans(polangle      )



# other settings for imshow
extent=( initial_xs[0], initial_xs[-1], initial_ys[0], initial_ys[-1])
interpolation = 'nearest'


plt.subplot(221)
plt.imshow(obs_hit_angle, interpolation=interpolation, extent=extent)
bar = plt.colorbar()
bar.formatter.set_useOffset(False)
bar.update_ticks()
plt.title('emitter angle')

plt.subplot(222)

plt.imshow(np.log10(redshift), interpolation=interpolation, extent=extent)
bar = plt.colorbar()
bar.formatter.set_useOffset(False)
bar.update_ticks()
plt.title('redshift')

plt.subplot(224)

plt.imshow(np.log10(bolflux), interpolation='lanczos', extent=extent)
bar = plt.colorbar()
bar.formatter.set_useOffset(False)
bar.update_ticks()
plt.title('bolometric flux')



plt.tight_layout()
plt.show()



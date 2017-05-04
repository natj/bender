import sys
sys.path.append('/Users/natj/projects/arcmancer/lib/')
import pyarcmancer as pyac
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'inferno'

#Setup star
R_in            = 12.0
m_in            = 1.4
freq_in         = 100.0
colat_in        = 90.0

#dist_in         = 100.0
#pts_per_axis_in = 20

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
planck_length_in_m = 1.616e-35
planck_mass_in_kg  = 2.176e-8
planck_time_in_s   = 5.3912e-44
planck_temp_in_K   = 1.417e32
planck_charge_in_C = 1.876e-18

# Operate in units where G = c = 1.
# Use units of solar masses
solar_mass_per_km = 0.6772
solar_mass_per_s  = 2.03e5


##################################################
#convert units to planck scale
#mass   = m_in * kg_per_Msun / planck_mass_in_kg
#rad    = R_in * 1000 / planck_length_in_m
#dist   = dist_in * 1000 / planck_length_in_m
#angvel = freq_in * 2.0*pi * planck_time_in_s
#colat  = colat_in * pi/180.0
#R_eq = rad


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
#conf.maximum_stepsize = 1.0
conf.maximum_steps     = 100000

conf.enforce_maximum_stepsize = False
conf.enforce_minimum_stepsize = True
conf.enforce_maximum_steps    = True
conf.store_only_endpoints     = True;


#Define metric of the spacetime
#metric = pyac.SchwarzschildMetric(R_eq, mass, angvel)
metric = pyac.AGMMetric(R_eq, 1.0, angvel, pyac.AGMMetric.MetricType.agm_standard)

#pyac.Log.set_console()
pyac.Log.set_file()


##################################################
# Construct image plane
#inclination = np.deg2rad(colat_in)
inclination = np.deg2rad(colat_in)
distance = mass * 1000

x_span = 6.0*R_eq
y_span = 6.0*R_eq
x_bins = 20
y_bins = 20

pixel_dx = 2*x_span / x_bins
pixel_dy = 2*y_span / y_bins
pixel_area = pixel_dx * pixel_dy


initial_xs = np.linspace(-x_span, x_span, x_bins)
initial_ys = np.linspace(-y_span, y_span, y_bins)


# construct coordinate tetrad
#def coordinate_tetrad(metric, position):
#    tetrad = [
#        np.array([1,0,0,0]), 
#        np.array([0,1,0,0]),
#        np.array([0,0,1,0]),
#        np.array([0,0,0,1]),
#    ]
#    return [ pyac.normalize(metric, position, v) for v in tetrad ]

## construct normalized coordinate vectors
#def coordinate_tetrad(metric, pos_bl):
#    tetrad = [
#        np.array([1.0,0,0,0]), 
#        np.array([0.0,1,0,0]),
#        np.array([0.0,0,1,0]),
#        np.array([0.0,0,0,1]),
#    ]
#    # gram-schmidt orthogonalize
#    for k in xrange(len(tetrad)):
#        uk = tetrad[k]
#        vk = tetrad[k]
#        #print "k = {}, u[{}] = {}".format(k, k, uk)
#        for j in xrange(k):
#            uj = tetrad[j]
#            #uk -= pyac.project_along(metric, pos_bl, tetrad[j], tetrad[k])
#            proj = metric.dot(pos_bl, vk, uj) / metric.dot(pos_bl, uj, uj) * uj
#            uk -= proj
#            #print "j = {}, proj = {}, u[{}] = {}".format(j, proj, k, uk)
#        tetrad[k] = uk
#    normalized_tetrad = [ pyac.normalize(metric, pos_bl, v) for v in tetrad ]
#    return normalized_tetrad

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


# generate an image plane of geodesics
def generate_image_plane(metric, distance, inclination, x_span, y_span, x_bins, y_bins):
    xs = initial_xs
    ys = initial_ys

    # avoid line of poles?
    xs += 0.1*pixel_dx

    plane = []
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            print "ix", ix, "x", x
            print "iy", iy, "y", y

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
                                np.dot(u_r     , vel_cart),
                                np.dot(u_theta , vel_cart),
                                np.dot(u_phi   , vel_cart)])

            # define vertical and horizontal
            vert = pyac.normalize(metric, x_sph, np.array([0, 0, -1.0, 0]))

            #vert_vel = pyac.project_along(metric, x_sph, vel_sph, vert)
            #vert -= vert_vel
            vert = pyac.normalize(metric, x_sph, vert)

            #horz = pyac.spatial_cross_product(
            #    metric, pyac.static_observer(metric, x_sph), vert, vel_sph)
            #horz = pyac.normalize(metric, x_sph, horz)
            horz = pyac.normalize(metric, x_sph, np.array([0, 0, 0, 1.0]))

            if 0:
                # test
                print "vert", vert, "horz", horz
                print "vert.horz", metric.dot(x_sph, vert, horz)
                print "vert.u", metric.dot(x_sph, vert, vel_sph)
                print "horz.u", metric.dot(x_sph, horz, vel_sph)


            # local polarization axes are: 
            # vertical = -u_theta
            # horizontal = u_phi
            #vert = np.array([0, 0, -1, 0])
            #horz = np.array([0, 0,  0, 1])

            geo = pyac.Geodesic(metric, x_sph, vel_sph, vert, horz, pyac.VectorType.null)
            plane.append((ix, iy, x, y, geo))


    return plane



##################################################
# AlGendy-Morsink metric specific stuff
def AGM_Ob(angvel, mass, rad):
    return 0.0
    #return angvel * np.sqrt(rad**3 / mass)

def AGM_o2(Ob, compactness):
    return 0.0
    
    #a0 = -0.788
    #a1 = 1.030
    #return Ob**2 * (a0 + a1*compactness)

def AGM_q(Ob, compactness):
    return 0.0
    # a2 = -0.11
    #return (Ob**2)* (a2/(compactness**2))

def AGM_beta(Ob, compactness):
    return 0.0
    #a1 = 0.4454
    #return a1 * (Ob**2)*compactness


def legendre_p2(x):
    return 0.5*(x**2 - 1.0)

def AGM_nu(mass, q, beta, x):
    r = x[2]
    th = x[3]
    nu0 = log(1.0 - mass/(2.0*r) / (1.0 + mass/(2.0*r)))
    return nu0 + (beta/3.0 - q*legendre_p2(np.cos(th))*mass/r)**3

def AGM_B(mass, beta, x):
    r = x[2]
    b0 = 1.0 - (mass/(2.0*r))**2
    return b0 + beta*(mass/r)**2


def AGM_isotropic2Sch(mass, q, beta, x):
    B = AGM_B(mass, beta, x)
    nu = AGM_nu(mass, q, beta, x)
    return B*exp(-nu)*x[2]


class NeutronStarSurface(pyac.Surface):

    #mass = 0.0
    #rad = 0.0
    #angvel = 0.0
    #compactness = 0.0
    #Omega_bar = 0.0
    #o2term = 0.0
    #q = 0.0
    #beta = 0.0

    @classmethod
    def from_arguments(cls, metric_, rad_, mass_, angvel_):
        ret = cls()
        ret.metric = metric_
        ret.mass = mass_
        ret.rad  = rad_
        ret.angvel = angvel_

        ret.compactness = ret.mass / ret.rad

        ret.Omega_bar = AGM_Ob(ret.angvel, ret.mass, ret.rad)
        ret.o2term = AGM_o2(ret.Omega_bar, ret.compactness)
        ret.q = AGM_q(ret.Omega_bar, ret.compactness)
        ret.beta = AGM_beta(ret.Omega_bar, ret.compactness)
        
        #ret.flattening = 1.0 - radius(0.0) / rad

        return ret

    def radius(self, colat):
        return self.rad * (1.0 + self.o2term * np.cos(colat)**2)

    def value(self, x):
        #r = x[2]

        print "self.mass=",self.mass
        print "self.q=",self.q
        print "self.beta=",self.beta

        rsurf = self.radius(x[3])
        r = AGM_isotropic2Sch(self.mass, self.q. self.beta, x)

        return rsurf - r

    def gradient(self, x): 
        #derivate of the radius function w.r.t. theta
        radius_th = self.rad * self.o2term * np.sin(2.0 * x[3])

        return np.array([0, 1, -radius_th, 0])


    def observer(self, metric, x):
        #return pyac.static_observer(metric, x)
        #return pyac.TangentVector(x, disk_fourvelocity(x, self.metric, self.mass, self.chi))

        #compute the four-velocity at the surface of a rotating
        #neutron star
        g = metric.value(x)

        #dphi / dt
        #= angular velocity as measured by observer at infinity
        #= (u^phi / u^t) evaluated at any radial distance
        dpdt = self.angvel
        
        #time component
        ut = 1.0/np.sqrt( g(0,0) + 2*dpdt * g(0, 3) + dpdt * dpdt * g(3,3) )

        #phi component
        uphi = dpdt * ut

        return pyac.TangentVector(x, np.array([ut, 0, 0, uphi]))



def compute_element(el, distance, metric, conf, surfaces):
    #print "computing pixel {},{}".format(el[0], el[1])
    el[4].compute(-2 * (np.tan(inclination)*x_span + distance), metric, conf, surfaces)
    #return el


# black-body specific intensity in natural units (h = 1)
def bb_intensity(nu, T):
    return 2 * nu**3 * ( np.exp( nu/T ) - 1 )**(-1)



##################################################
#x0 = np.array([0, 1, 0, 0])
#v0 = np.array([0, 1, 0, 0])
#geo = pyac.Geodesic(metric, x0, v0, pyac.VectorType.null)
#print metric.sqnorm(geo.get_points()[0].point), -1.0
#
###################################################
#geo = pyac.Geodesics(metric, 
#
#geo.compute(10, metric, conf)
#print geo.get_points()[-1].point.x[1], x0[1] + 10*v0[1]
#print geo.get_points()



imgplane = generate_image_plane(metric, distance, inclination, x_span,
                                y_span, x_bins, y_bins)


print R_eq
print mass
print angvel

#ns_surface = NeutronStarSurface.from_arguments(metric, R_eq, mass, angvel)

#ns_surface = pyac.AGMSurface(R_eq, mass, angvel, pyac.AGMSurface.SurfaceType.spherical)
ns_surface = pyac.AGMSurface(R_eq, mass, angvel, pyac.AGMSurface.SurfaceType.oblate)

surfaces = [ ns_surface ]

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



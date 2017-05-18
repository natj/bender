import sys
sys.path.append('/Users/natj/projects/arcmancer/lib/')
import pyarcmancer as pyac

import numpy as np
import scipy.interpolate as interp


##################################################
# Auxiliary functions


#modulo 2pi
def mod2pi(x):
    while (x > 2.0*np.pi):
        x -= 2.0*np.pi
    while (x < 0.0):
        x += 2.0*np.pi
    return x

# Transform cartesian xy to our special polar angle
def calc_chi(x,y):
    return mod2pi(np.pi/2.0 - np.arctan2(y,x) )

# (x,y) to (r, chi)
def xy2pol(x, y):
    return np.hypot(x, y), calc_chi(x, y)





class Imgplane:


    incl     = 90.0
    distance = 150.0
    rmax     = 10.0
    time0    = 0.0

    verbose  = 0

    def __init__(self, conf, metric, surfaces):
        self.conf  = conf
        self.metric = metric
        self.surfaces = surfaces


    # in the limit m -> 0, BL coordinates go to oblate spheroidal minkowski. These
    # go to minkowski for r -> inf or a ->0
    def boyer_lindquist_position(self, x_cart):
        #print "Transforming cartesian position {} to Boyer-Lindquist".format(x_cart)
        x, y, z = x_cart
        r = np.linalg.norm(x_cart)
        theta = np.arccos(z/r)
        phi = np.arctan2(y, x)
        return np.array([0, r, theta, phi])


    # construct local spherical axes in cartesian coordinates
    def local_spherical_axes(self, pos_cart):
        #print "Computing cartesian components of local spherical axes at {}".format(pos_cart)
        u_r = pos_cart / np.linalg.norm(pos_cart)
        u_phi = np.cross(np.array([0,0,1]), u_r)
        u_phi /= np.linalg.norm(u_phi)
        u_theta = -np.cross(u_r, u_phi)
        u_theta /= np.linalg.norm(u_theta)
        #print "result {} {} {}".format(u_r, u_theta, u_phi)
        return [u_r, u_theta, u_phi]

    #Build geodesic path from xy coordinates in image plane
    def xy2geo(self, x, y):

        # get coordinates for position of image plane point
        normal = np.array([np.sin(self.incl), 
                           0.0, 
                           np.cos(self.incl)]
                         )

        x_cart = np.array([-y*np.cos(self.incl), 
                           x, 
                           y*np.sin(self.incl)])

        #shift to distance
        x_cart += self.distance * normal
        x_sph = self.boyer_lindquist_position(x_cart)


        # get velocity by projecting to local B-L axes
        # get axes in _cartesian_ coordinates
        u_r, u_theta, u_phi = self.local_spherical_axes(x_cart)
        vel_cart = normal
        vel_sph = np.array([0, 
                            np.dot(u_r     , vel_cart) ,
                            np.dot(u_theta , vel_cart) / x_sph[1],
                            np.dot(u_phi   , vel_cart) / (x_sph[1] * np.sin(x_sph[2]))])


        # define vertical and horizontal
        vert = pyac.normalize(self.metric, x_sph, np.array([0, 0, -1.0, 0]))

        vert_vel = pyac.project_along(self.metric, x_sph, vel_sph, vert)
        vert -= vert_vel
        vert = pyac.normalize(self.metric, x_sph, vert)

        horz = pyac.spatial_cross_product(
            self.metric, pyac.static_observer(self.metric, x_sph), vert, vel_sph)
        horz = pyac.normalize(self.metric, x_sph, horz)
        horz = pyac.normalize(self.metric, x_sph, np.array([0, 0, 0, 1.0]))


        geo = pyac.Geodesic(self.metric, x_sph, vel_sph, vert, horz, pyac.VectorType.null)

        return geo




    def polar2geo(self, rad, chi):
        x = rad * np.sin(chi)
        y = rad * np.cos(chi)

        return self.xy2geo(x, y)


    def compute_element(self, geo):
        geo.compute(-(self.distance + 50.0), self.metric, self.conf, self.surfaces)


    #Rotate in chi angle and find star boundaries
    def find_boundaries(self, 
                        Nedge=10,
                        reltol = 1.0e-3,
                        max_iterations = 20
                        ):
        if self.verbose > 0:
            print "Finding edge boundaries for the star..."
            print "  # angles       {}".format(Nedge)
            print "  max iterations {}".format(max_iterations)
            print "  relative tol   {}".format(reltol)


        chis = np.linspace(0.0, 1.0, Nedge)*2.0*np.pi + 0.001
        rlims = np.zeros(Nedge)

        rmin = 0.0
        rmax = 12.0


        for i, chii in enumerate(chis):
            geos = []

            rmini = rmin
            rmaxi = rmax
            rmid = 0.0

            N = 0

            relerr = 1.0
            rmid_old = 100.0

            while (N < max_iterations) and (relerr > reltol):
                rmid = (rmini + rmaxi)/2.0

                geos.append(self.polar2geo( rmid, chii) )
                self.compute_element(geos[N])

                hit = geos[N].front_termination().hit_surface

                if hit:
                    rmini = rmid
                else:
                    rmaxi = rmid
                
                relerr = np.abs(rmid - rmid_old)/rmid
                rmid_old = rmid

                N += 1
                if self.verbose > 1:
                    print "Iterating edge at {} after {} tries for angle={} ({})".format(rmid, N, chii, relerr)

            rlims[i] = rmid


        self.rmax = np.max(rlims)*1.01
        if self.verbose > 1:
            print "Maximum edge {}".format(self.rmax)

        #Build edge location interpolator
        self.edge = interp.InterpolatedUnivariateSpline(chis, rlims)

         

    ##################################################
    # Creates internal polar grid for the interpolation
    def generate_internal_grid(self, 
                            Nrad = 30, 
                            Nchi = 20,
                            use_flat_chi = False
                            ):

        if self.verbose > 0:
            print "Generating internal polar grid..."

        self.Nrad = Nrad
        self.Nchi = Nchi

        
        #Build non-equidistant angle grid 
        # we specifically add points near the chi = 0
        # and ease the interpolation with two egde points 
        dchi_edge = 0.001
        chimin    = 0.0 - dchi_edge
        chimax    = 2.0*np.pi + dchi_edge
        
        chi_diffs     = 0.8 + np.sin( np.linspace(0.0, 2.0*np.pi, Nchi-3) )**2
        chi_diffs     = np.insert(chi_diffs, 0, 0.0)
        self.chi_grid = chimin + (chimax - chimin) * np.cumsum(chi_diffs)/np.sum(chi_diffs)
        self.chi_grid = np.insert(self.chi_grid, 0, self.chi_grid[0] - dchi_edge)
        self.chi_grid = np.append(self.chi_grid, self.chi_grid[-1] + dchi_edge)
    

        #If true, we fall back to linear angle spacing
        if use_flat_chi:
            self.chi_grid = np.linspace(0.0, 2.0*np.pi, Nchi)
        
        
        # Build non-equidistant radius grid
        # The grid is Gauss-Laguerre weighted so that
        # more points are near the boundary where the
        # surface curves more
        #
        # Actual computations are black magic so only
        # modify these if you know what you are doing
        rad_diffs     = 1.0 / np.exp( np.linspace(1.2, 2.0, Nrad-1)**2)
        self.rad_grid = self.rmax * np.cumsum(rad_diffs)/np.sum(rad_diffs)
        self.rad_grid = np.insert(self.rad_grid, 0, 0.001)
        
        
        self.grid = np.empty((Nrad,Nchi), dtype=np.object)
        
        for i, chi in enumerate(self.chi_grid):
            if self.verbose > 0:
                if i % 10 == 0:
                    print "{} % done".format(float(i)/len(self.chi_grid) * 100)

            for j, rad in enumerate(self.rad_grid):
                if self.verbose > 2:
                    print "  tracing geodesic at chi={} and r={}".format(chi, rad) 

                #Trace geodesic from image plane to star
                geo = self.polar2geo(rad, chi)
                self.compute_element(geo)

                self.grid[j,i] = geo
        


    #Reduce geodesic into physical quantities
    def dissect_geos(self):
    
        if self.verbose > 0:
            print "Dissecting geodesic paths to observed quantities..."
    
        Nrad = self.Nrad
        Nchi = self.Nchi
    
        self.Reds   = np.zeros((Nrad, Nchi))
        self.Cosas  = np.zeros((Nrad, Nchi))
        self.Times  = np.zeros((Nrad, Nchi))
        self.Thetas = np.zeros((Nrad, Nchi))


        self.Phis_sin   = np.zeros((Nrad, Nchi))
        self.Phis_cos   = np.zeros((Nrad, Nchi))
    
        
        for i, chi in enumerate(self.chi_grid):
            if self.verbose > 1:
                print "{} % done".format(float(i)/len(chi_grid) * 100)

            for j, rad in enumerate(self.rad_grid):
    
                geo = self.grid[j,i]
    
                hit_pt = geo.get_points()[0]
                obs_pt = geo.get_points()[-1]
    
                #skip if we did not hit
                hit = geo.front_termination().hit_surface

                if not(hit):
                    continue


                #coordinates 
                #surface_point.point.x
                t  = self.Times[j,i]  = hit_pt.point.x[0]
                th = self.Thetas[j,i] = hit_pt.point.x[2]

                ps  = self.Phis_sin[j,i] = np.sin(hit_pt.point.x[3])
                pc  = self.Phis_cos[j,i] = np.cos(hit_pt.point.x[3])
                
    
                #Redshift
                g = self.Reds[j,i] = \
                    self.metric.dot(obs_pt.point, pyac.static_observer(self.metric, obs_pt.x())) / \
                        self.metric.dot(hit_pt.point, self.surfaces[0].observer(self.metric, hit_pt.x()))
    
                #hit angle
                cosa = self.Cosas[j,i] = \
                    geo.front_termination().observer_hit_angle
    


        #Now compute spline interpolators 
        if self.verbose > 0:
            print "Building spline coefficients..."

        kx = ky = 2
        s = 0


        #Shift time XXX
        time_, phi_, theta_, cosa_, reds_ = self.get_exact_pixel(0.0, 0.0)
        self.time0 = time_
        self.Times -= time_


        self.intp_Times = interp.RectBivariateSpline(self.rad_grid, self.chi_grid, self.Times, kx=kx, ky=ky, s=s)

        self.intp_Thetas = interp.RectBivariateSpline(self.rad_grid, self.chi_grid, self.Thetas, kx=kx, ky=ky, s=s)


        self.intp_Phis_sin = interp.RectBivariateSpline(self.rad_grid, self.chi_grid, self.Phis_sin, kx=kx, ky=ky, s=s)
        self.intp_Phis_cos = interp.RectBivariateSpline(self.rad_grid, self.chi_grid, self.Phis_cos, kx=kx, ky=ky, s=s)


        self.intp_Cosas = interp.RectBivariateSpline(self.rad_grid, self.chi_grid, self.Cosas, kx=kx, ky=ky, s=s)

        self.intp_Reds = interp.RectBivariateSpline(self.rad_grid, self.chi_grid, self.Reds, kx=kx, ky=ky, s=s)



    
    #TODO define return class
    def get_pixel(self, x, y):
        rad, chi = xy2pol(x, y)

        time = 0.0
        phi  = 0.0
        theta = 0.0
        cosa = 0.0
        reds = 0.0


        if rad <= self.edge(chi):
            time  = self.intp_Times.ev(rad, chi)

            phis  = self.intp_Phis_sin.ev(rad, chi)
            phic  = self.intp_Phis_cos.ev(rad, chi)
            phi   = np.arctan2(phis, phic)

            theta = self.intp_Thetas.ev(rad, chi)
            cosa  = self.intp_Cosas.ev(rad, chi)
            reds  = self.intp_Reds.ev(rad, chi)

            #pix = pixel([time, theta, phi], reds, cosa)

        return time, phi, theta, cosa, reds
    
    def get_exact_pixel(self, x, y):

        geo = self.xy2geo(x, y)
        self.compute_element(geo)

        hit = geo.front_termination().hit_surface
        if not(hit):
            return 0.0, 0.0, 0.0, 0.0, 0.0

        hit_pt = geo.get_points()[0]
        obs_pt = geo.get_points()[-1]
    
        #coordinates 
        time  = hit_pt.point.x[0] - self.time0
        theta = hit_pt.point.x[2]
        phi   = hit_pt.point.x[3]

        #Redshift
        reds = self.metric.dot(obs_pt.point, pyac.static_observer(self.metric, obs_pt.x())) / \
            self.metric.dot(hit_pt.point, self.surfaces[0].observer(self.metric, hit_pt.x()))
    
        #hit angle
        cosa = geo.front_termination().observer_hit_angle
    
        return time, phi, theta, cosa, reds



# Pixel class that connects the image plane coordinates to star coordinates
class pixel:

    def __init__(self, 
                 x, 
                 redshift, 
                 hit_angle):

        #self.x = np.array([0.0, 0.0, 0.0])
        self.x = x

        self.redshift  = redshift
        self.hit_angle = hit_angle






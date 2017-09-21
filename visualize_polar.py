import numpy as np
import matplotlib as mpl
from pylab import *
from matplotlib import cm



##################################################
# plot values on image plane
def trans(mat):
    #return np.flipud(mat.T) #trans
    return np.flipud(mat).T  #detrans

def detrans(mat):
    return mat  #detrans


#mask all 0.0 elements and self.origin
def mask(mat):
    return np.ma.masked_where(mat == 0, mat) 



#build up a chess board layering from phi and theta coordinates
def chess_layer(phis, thetas, star_rot):

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

            xd = np.int(60.0*(phi + star_rot)/(2*pi))
            yd = np.int(60.0*theta/(2*pi))
            
            if (xd & 1) and (yd & 1):
                mat[i,j] = black
            elif (xd & 0) and (yd & 0):
                mat[i,j] = black
            else:
                mat[i,j] = white

    return mat



##################################################
#Visualization class for neutron star figures
class Visualize:

    #gridspec size
    nx = 3
    ny = 6

    #image plane width & height
    x_span = 10.0
    y_span = 10.0
    
    x_bins = 500
    y_bins = 500

    px_bins = 200
    py_bins = 200

    #Interpolation type for the imshow 
    interpolation = 'nearest'

    #Compactness for color scalings
    compactness = 1.0

    #flag for not re-plotting everything
    initialized = False

    #imshow parameters
    origin='lower'

    def __init__(self):

        self.gs = GridSpec(self.ny, self.nx)
        self.gs.update(hspace = 0.5)
        self.axs = np.empty(( self.nx*self.ny ), dtype=np.object)
        
        #Big neutron star visualization panel
        self.axs[0] = subplot( self.gs[0:2, 0:2] )
        self.axs[0].axis('off')


        #create line object ready for spot bounding box
        #self.line0, = self.axs[0].plot([0,0,0,0,0],[0,0,0,0,0],"k-")


        #Other minor figures
        self.axs[1] = subplot( self.gs[2, 0] )
        self.axs[2] = subplot( self.gs[2, 1] )

        self.axs[3] = subplot( self.gs[0, 2] )
        self.axs[4] = subplot( self.gs[1, 2] )
        self.axs[5] = subplot( self.gs[2, 2] )


        for i in range(1,6):
            self.axs[i].minorticks_on()

        
        self.axs[1].set_title(r'emitter angle $\alpha$')
        self.axs[2].set_title(r'redshift')
        self.axs[3].set_title(r'$\phi$')
        self.axs[4].set_title(r'$\theta$')
        self.axs[5].set_title(r'time')
          

        self.pixel_dx = 2*self.x_span / self.x_bins
        self.pixel_dy = 2*self.y_span / self.y_bins
        self.pixel_area = self.pixel_dx * self.pixel_dy

        self.xs = np.linspace(-self.x_span, self.x_span, self.x_bins)
        self.ys = np.linspace(-self.y_span, self.y_span, self.y_bins)


        #build interpolated array
        self.spotarea      = np.zeros((self.x_bins, self.y_bins))
        self.redshift      = np.zeros((self.x_bins, self.y_bins))
        self.obs_hit_angle = np.zeros((self.x_bins, self.y_bins))
        self.times         = np.zeros((self.x_bins, self.y_bins))
        self.thetas        = np.zeros((self.x_bins, self.y_bins))
        self.phis          = np.zeros((self.x_bins, self.y_bins))


        # other settings for imshow
        self.extent=( self.xs[0], self.xs[-1], self.ys[0], self.ys[-1] )


        #polar grid
        self.axs[7] = subplot( self.gs[4:6,:] )
        self.axs[7].minorticks_on()
        self.axs[7].set_xlabel(r'$\chi$ $[R_{max}]$ ')
        self.axs[7].set_ylabel(r'$r$ $[\pi]$ ')

        self.pxs = np.linspace(0.0, 2.0, self.px_bins)
        self.pys = np.linspace(0.0, 1.0, self.py_bins)

        self.pextent=( self.pxs[0], self.pxs[-1], self.pys[0], self.pys[-1] )

        #build interpolated array
        self.polar_spotarea      = np.zeros((self.px_bins, self.py_bins))
        self.polar_redshift      = np.zeros((self.px_bins, self.py_bins))
        self.polar_obs_hit_angle = np.zeros((self.px_bins, self.py_bins))
        self.polar_times         = np.zeros((self.px_bins, self.py_bins))
        self.polar_thetas        = np.zeros((self.px_bins, self.py_bins))
        self.polar_phis          = np.zeros((self.px_bins, self.py_bins))



    #Get observables from img
    def dissect(self, img):

        for i, xi in enumerate(self.xs):
            for j, yi in enumerate(self.ys):

                time, phi, theta, cosa, reds = img.get_pixel(xi, yi)
                #time, phi, theta, cosa, reds = img.get_exact_pixel(xi, yi)
        
                self.redshift[i,j]      = reds
                self.obs_hit_angle[i,j] = cosa
                self.times[i,j]         = time
                self.thetas[i,j]        = theta
                self.phis[i,j]          = phi



    #Get observables from img
    def polar_dissect(self, img):

        for i, c in enumerate(self.pxs):
            for j, r in enumerate(self.pys):

                #rad = r*img.rmax
                chi = c*np.pi
                rad = img.radius_stretch(r, chi) #* img.rmax

                time, phi, theta, cosa, reds  = img.get_poxel(rad, chi)
                self.polar_redshift[i,j]      = reds
                self.polar_obs_hit_angle[i,j] = cosa
                self.polar_times[i,j]         = time
                self.polar_thetas[i,j]        = theta
                self.polar_phis[i,j]          = phi





    #Separate dissection for patterns on the surface
    def polar_dissect_spot(self, img, spot):
        for i, c in enumerate(self.pxs):
            for j, r in enumerate(self.pys):

                #rad = r*img.rmax
                #chi = c*np.pi

                time  = self.polar_times[i,j]
                phi   = self.polar_phis[i,j]
                theta = self.polar_thetas[i,j]
        
                hits_spot = spot.hit([time, theta, phi])

                #Change state if we found the spot
                #if hits_spot:

                #record hit pixels
                self.polar_spotarea[i,j] = 1.0 if hits_spot else 0.0




    #Separate dissection for patterns on the surface
    def dissect_spot(self, img, spot):

        self.frame_x1 = self.xs[-1]
        self.frame_x2 = self.xs[0]

        self.frame_y1 = self.ys[-1]
        self.frame_y2 = self.ys[0]
        
        self.located_spot = False

        for i, xi in enumerate(self.xs):
            for j, yi in enumerate(self.ys):

                time  = self.times[i,j]
                phi   = self.phis[i,j]
                theta = self.thetas[i,j]

                if self.redshift[i,j] == 0.0:
                    continue
        
                hits_spot = spot.hit([time, theta, phi])

                #Change state if we found the spot
                if hits_spot:
                    self.located_spot = True

                    #record bounding boxes
                    self.frame_y2 = yi if self.frame_y2 < yi else self.frame_y2 #top max
                    self.frame_y1 = yi if self.frame_y1 > yi else self.frame_y1 #bot min

                    self.frame_x1 = xi if self.frame_x1 > xi else self.frame_x1 #left min
                    self.frame_x2 = xi if self.frame_x2 < xi else self.frame_x2 #right max


                #record hit pixels
                self.spotarea[i,j] = 1.0 if hits_spot else 0.0



    # Locate spot boundaries by projecting spot array to x/y axis 
    # then look for edges of histogram
    def improve_on_spot_boundaries(self):

        yproj = np.sum(self.spotarea, 0)
        xproj = np.sum(self.spotarea, 1)

        xmi = np.argmax(xproj)
        ymi = np.argmax(yproj)
        xm = self.xs[xmi]
        ym = self.ys[ymi]

        ymin = 0.0
        ymax = 0.0

        i = ymi
        while i > 0:
            if yproj[i] == 0.0:
                break
            ymin = self.ys[i]
            i -= 1

        i = ymi
        while i < len(self.ys):
            if yproj[i] == 0.0:
                break
            ymax = self.ys[i]
            i += 1

        self.frame_y1 = ymin
        self.frame_y2 = ymax

        xmin = 0.0
        xmax = 0.0

        i = xmi
        while i > 0:
            if xproj[i] == 0.0:
                break
            xmin = self.xs[i]
            i -= 1

        i = xmi
        while i < len(self.xs):
            if xproj[i] == 0.0:
                break
            xmax = self.xs[i]
            i += 1

        self.frame_x1 = xmin
        self.frame_x2 = xmax


    ################################################## 
    #Actual visualizations

    #Star plot
    def star_plot(self, star_rotation):

        #Compute chess pattern
        chess = chess_layer(self.phis, self.thetas, star_rotation)
        chess    = trans(mask(chess))

        redshift = trans(mask(self.redshift))
        spotarea = trans(mask(self.spotarea))

        self.axs[0].clear()
        self.axs[0].axis('off')

        self.axs[0].imshow(
                chess,
                interpolation=self.interpolation, 
                extent=self.extent, 
                cmap=cm.get_cmap('Greys'), 
                vmin=0.8, 
                vmax=2.0, 
                alpha=0.8, 
                origin=self.origin
                )


        self.axs[0].imshow(
                redshift,
                interpolation=self.interpolation, 
                origin=self.origin,
                extent=self.extent,
                cmap=cm.get_cmap('coolwarm_r'), 
                vmin=0.9*self.compactness, 
                vmax=1.1*self.compactness, 
                alpha=0.95
                )
    
        levels = np.linspace(0.8*self.compactness, 1.2*self.compactness, 20)

        self.axs[0].contour(
                redshift,
                levels, 
                hold='on', 
                colors='w',
                origin=self.origin,
                extent=self.extent, 
                vmin=0.8*self.compactness, 
                vmax=1.2*self.compactness
                )
    
        self.axs[0].imshow(
                spotarea,
                interpolation=self.interpolation, 
                extent=self.extent, 
                cmap=cm.get_cmap('inferno'), 
                alpha=0.7,
                origin=self.origin
                )
    

    def polar_plot(self, star_rotation):

        #Compute chess pattern
        chess = chess_layer(self.polar_phis, self.polar_thetas, star_rotation)
        chess    = trans(mask(chess))

        redshift = trans(mask(self.polar_redshift))
        spotarea = trans(mask(self.polar_spotarea))

        self.axs[7].clear()
        #self.axs[6].axis('off')

        self.axs[7].imshow(
                chess,
                interpolation=self.interpolation, 
                extent=self.pextent, 
                cmap=cm.get_cmap('Greys'), 
                vmin=0.8, 
                vmax=2.0, 
                alpha=0.8, 
                aspect='auto',
                origin=self.origin
                )


        self.axs[7].imshow(
                redshift,
                interpolation=self.interpolation, 
                origin=self.origin,
                extent=self.pextent,
                cmap=cm.get_cmap('coolwarm_r'), 
                vmin=0.9*self.compactness, 
                vmax=1.1*self.compactness, 
                aspect='auto',
                alpha=0.95
                )
    
        levels = np.linspace(0.8*self.compactness, 1.2*self.compactness, 20)

        self.axs[7].contour(
                redshift,
                levels, 
                hold='on', 
                colors='w',
                origin=self.origin,
                extent=self.pextent, 
                aspect='auto',
                vmin=0.8*self.compactness, 
                vmax=1.2*self.compactness
                )
    
        self.axs[7].imshow(
                spotarea,
                interpolation=self.interpolation, 
                extent=self.pextent, 
                cmap=cm.get_cmap('inferno'), 
                alpha=0.7,
                aspect='auto',
                origin=self.origin
                )
    


    #Plot img class & spot
    def star(self, img, spot):

        self.dissect(img)
        self.dissect_spot(img, spot)

        phirot = -spot.star_time * spot.angvel
        self.star_plot(phirot)


    def polar(self, img, spot):

        self.polar_dissect(img)
        self.polar_dissect_spot(img, spot)

        phirot = -spot.star_time * spot.angvel
        self.polar_plot(phirot)



    def plot(self, img):

        self.dissect(img)

        redshift      = trans(mask(self.redshift))
        obs_hit_angle = trans(mask(self.obs_hit_angle))
        times         = trans(mask(self.times))
        thetas        = trans(mask(self.thetas))
        phis          = trans(mask(self.phis))


        #observer hit angle
        cax = self.axs[1].imshow(
                obs_hit_angle, 
                interpolation=self.interpolation, 
                extent=self.extent,
                origin=self.origin
                )
        #colorbar(cax)

        #redshift
        cax = self.axs[2].imshow(
                redshift, 
                interpolation=self.interpolation, 
                origin=self.origin,
                extent=self.extent,
                cmap=cm.get_cmap('coolwarm_r'), 
                vmin=0.85*self.compactness, 
                vmax=1.15*self.compactness
                )
        
        levels = np.linspace(0.8*self.compactness, 1.2*self.compactness, 20)
        self.axs[2].contour(
                redshift, 
                levels, 
                hold='on', 
                colors='w',
                origin=self.origin,
                extent=self.extent, 
                vmin=0.8*self.compactness, 
                vmax=1.2*self.compactness
                )
        #colorbar(cax)


        #phi angle
        cax = self.axs[3].imshow(
                phis, 
                interpolation=self.interpolation, 
                extent=self.extent,
                origin=self.origin
                )
        #colorbar(cax)

        #theta angle
        cax = self.axs[4].imshow(
                thetas, 
                interpolation=self.interpolation, 
                extent=self.extent,
                origin=self.origin
                )
        #colorbar(cax)

        cax = self.axs[5].imshow(
                times, 
                interpolation=self.interpolation, 
                extent=self.extent,
                origin=self.origin
                )
        #levels = np.linspace(0.8*self.compactness, 1.2*self.compactness, 20)
        self.axs[5].contour(
                times, 
                10, 
                hold='on', 
                colors='w',
                origin=self.origin,
                extent=self.extent, 
                #vmin=0.8*self.compactness, 
                #vmax=1.2*self.compactness
                )




    #Show cartesian bounding box surrounding the box
    def spot_bounding_box(self):

        #expand image a bit
        frame_expansion_x = 1.0*np.abs(self.xs[1] - self.xs[0])
        frame_expansion_y = 1.0*np.abs(self.ys[1] - self.ys[0])

        self.frame_x1 -= frame_expansion_x
        self.frame_x2 += frame_expansion_x
        self.frame_y1 -= frame_expansion_y
        self.frame_y2 += frame_expansion_y

        curvex = [-self.frame_x1,
                  -self.frame_x2,
                  -self.frame_x2,
                  -self.frame_x1,
                  -self.frame_x1]
        curvey = [self.frame_y1,
                  self.frame_y1,
                  self.frame_y2,
                  self.frame_y2,
                  self.frame_y1]

        if self.located_spot:
            self.axs[0].plot(curvex, curvey, "k-")
        

        return [self.frame_x1, self.frame_y1, self.frame_x2, self.frame_y2]
        

    




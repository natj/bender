import numpy as np
import matplotlib as mpl
from pylab import *
from matplotlib import cm






##################################################
# plot values on image plane
def trans(mat):
    return np.flipud(mat.T)
#return mat
def detrans(mat):
    return np.flipud(mat).T

#mask all 0.0 elements and transpose
def mask(mat):
    return np.ma.masked_where(mat == 0, mat) 



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



##################################################
#Visualization class for neutron star figures
class Visualize:

    #gridspec size
    nx = 3
    ny = 3

    #image plane width & height
    x_span = 10.0
    y_span = 10.0
    
    x_bins = 200
    y_bins = 200

    #Interpolation type for the imshow 
    interpolation = 'nearest'

    #Compactness for color scalings
    compactness = 1.0

    #flag for not re-plotting everything
    initialized = False

    def __init__(self):

        self.gs = GridSpec(self.ny, self.nx)
        self.axs = np.empty(( self.nx*self.ny ), dtype=np.object)
        
        #Big neutron star visualization panel
        self.axs[0] = subplot( self.gs[0:2, 0:2] )
        self.axs[0].axis('off')


        #Other minor figures
        self.axs[1] = subplot( self.gs[2, 0] )
        self.axs[2] = subplot( self.gs[2, 1] )

        self.axs[3] = subplot( self.gs[0, 2] )
        self.axs[4] = subplot( self.gs[1, 2] )

        for i in range(1,5):
            self.axs[i].minorticks_on()

        
        self.axs[1].set_title(r'emitter angle $\alpha$')
        self.axs[2].set_title(r'redshift')
        self.axs[3].set_title(r'$\phi$')
        self.axs[4].set_title(r'$\theta$')
          


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
        self.extent=( self.xs[0], self.xs[-1], self.ys[0], self.xs[-1] )


    #Get observables from img
    def dissect(self, img):

        for i, xi in enumerate(self.xs):
            for j, yi in enumerate(self.ys):
                time, phi, theta, cosa, reds = img.get_pixel(xi, yi)
        
                self.redshift[i,j]      = reds
                self.obs_hit_angle[i,j] = cosa
                self.times[i,j]         = time
                self.thetas[i,j]        = theta
                self.phis[i,j]          = phi

        self.redshift      = trans(self.redshift)
        self.obs_hit_angle = trans(self.obs_hit_angle)
        self.times         = trans(self.times)
        self.thetas        = trans(self.thetas)
        self.phis          = trans(self.phis)



    #Separate dissection for patterns on the surface
    def dissect_spot(self, img, spot):

        for i, xi in enumerate(self.xs):
            for j, yi in enumerate(self.ys):

                time  = self.times[i,j]
                phi   = self.phis[i,j]
                theta = self.thetas[i,j]
        
                hits_spot = spot.hit([time, theta, phi])
        
                self.spotarea[i,j] = 1.0 if hits_spot else 0.0



    ################################################## 
    #Actual visualizations

    #Star plot
    def star(self):

        #Compute chess pattern
        chess = chess_layer(self.phis, self.thetas)

        self.axs[0].imshow(
                mask(chess), 
                interpolation=self.interpolation, 
                extent=self.extent, 
                cmap=cm.get_cmap('Greys'), 
                vmin=0.8, 
                vmax=2.0, 
                alpha=0.6, 
                origin='lower')


        self.axs[0].imshow(
                mask(self.redshift),
                interpolation=self.interpolation, 
                origin='lower', 
                extent=self.extent,
                cmap=cm.get_cmap('coolwarm_r'), 
                vmin=0.9*self.compactness, 
                vmax=1.1*self.compactness, 
                alpha=0.95
                )
    
        levels = np.linspace(0.9*self.compactness, 1.2*self.compactness, 20)

        self.axs[0].contour(
                mask(self.redshift),
                levels, 
                hold='on', 
                colors='w',
                origin='lower', 
                extent=self.extent, 
                vmin=0.8*self.compactness, 
                vmax=1.2*self.compactness
                )
    
        self.axs[0].imshow(
                mask(self.spotarea),
                interpolation=self.interpolation, 
                extent=self.extent, 
                cmap=cm.get_cmap('inferno'), 
                alpha=0.7
                )
    


    #Plot img class & spot
    def plot_star(self, img, spot):

        self.dissect(img)
        self.dissect_spot(img, spot)

        self.star()
        pause(0.001)


    def plot(self, img):

        self.dissect(img)

        #observer hit angle
        cax = self.axs[1].imshow(
                mask(self.obs_hit_angle), 
                interpolation=self.interpolation, 
                extent=self.extent,
                origin='lower'
                )
        #colorbar(cax)

        #redshift
        cax = self.axs[2].imshow(
                mask(self.redshift), 
                interpolation=self.interpolation, 
                origin='lower', 
                extent=self.extent,
                cmap=cm.get_cmap('coolwarm_r'), 
                vmin=0.85*self.compactness, 
                vmax=1.15*self.compactness
                )
        
        levels = np.linspace(0.9*self.compactness, 1.1*self.compactness, 20)
        self.axs[2].contour(
                mask(self.redshift), 
                levels, 
                hold='on', 
                colors='w',
                origin='lower', 
                extent=self.extent, 
                vmin=0.85*self.compactness, 
                vmax=1.15*self.compactness
                )
        #colorbar(cax)


        #phi angle
        cax = self.axs[3].imshow(
                mask(self.phis), 
                interpolation=self.interpolation, 
                extent=self.extent,
                origin='lower'
                )
        #colorbar(cax)

        #theta angle
        cax = self.axs[4].imshow(
                mask(self.thetas), 
                interpolation=self.interpolation, 
                extent=self.extent,
                origin='lower'
                )
        #colorbar(cax)


        pause(0.001)





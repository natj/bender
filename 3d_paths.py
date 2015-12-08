import numpy as np
import math
from pylab import *
import csv
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


from os import listdir
from os.path import isfile, join



def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]; x_mean = np.mean(x_limits)
    y_range = y_limits[1] - y_limits[0]; y_mean = np.mean(y_limits)
    z_range = z_limits[1] - z_limits[0]; z_mean = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mean - plot_radius, x_mean + plot_radius])
    ax.set_ylim3d([y_mean - plot_radius, y_mean + plot_radius])
    ax.set_zlim3d([z_mean - plot_radius, z_mean + plot_radius])

    

def read_path(fname):
    da = np.genfromtxt(fname, delimiter=',')
    return da[:,0], da[:,1], da[:,2], da[:,3], da[:,4]


def spher_2_cart(rad, theta, phi, rcutmin=0.0, rcutmax=20.0):

    xs = []
    ys = []
    zs = []
    
    N = len(rad)
    for i in range(N):
        if rad[i] != 0.0:
            if rcutmin < (1/rad[i]) < rcutmax:
                rad1 = rad[i]
                the1 = theta[i]
                phi1 = phi[i]
                #print rad1, the1, phi1
            
                x = (1/rad1)*cos(phi1)*sin(the1)
                y = (1/rad1)*sin(phi1)*sin(the1)
                z = (1/rad1)*cos(the1)

                #print x, y, z
                xs.append(x)
                ys.append(y)
                zs.append(z)
                
    #return xs,ys,zs
    return np.asarray(xs),np.asarray(ys),np.asarray(zs)

            
#fdir = "ppath_x2_yslice/"
#fdir = "ppath/"
fdir = "ppath_inclpi4/"

#####################
fig = figure()
rc('font', family='serif')
rc('xtick', labelsize='x-small')
rc('ytick', labelsize='x-small')

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")



#plot initial sphere-wireframe
R = 3.4
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
#x = R * np.outer(np.cos(u), np.sin(v))
#y = R * np.outer(np.sin(u), np.sin(v))
#z = R * np.outer(np.ones(np.size(u)), np.cos(v))

#oblate spheroid
ecc = 0.648
eta = 0.5*np.log((2-ecc**2)/ecc**2)
x = R * np.cosh(eta)*np.outer(np.cos(u), np.sin(v))
y = R * np.cosh(eta)*np.outer(np.sin(u), np.sin(v))
z = R * np.sinh(eta)*np.outer(np.ones(np.size(u)), np.cos(v))


ax.plot_surface(x, y, z,
                rstride=4, cstride=4,
                alpha = 0.05,
                color='b')


#plot photons
onlyfiles = [f for f in listdir(fdir) if isfile(join(fdir, f))]


for i in range(len(onlyfiles)):
    pfile = onlyfiles[i]

    if 'p_' in pfile:
        fname = fdir+pfile

        #extract x and y from dirname
        xx = float(pfile.split('_')[1])
        yy = float((pfile.split('_')[2]).split('.csv')[0])
                
        print fname, xx, yy

        rad, theta, phi, err, lvl = read_path(fname)
        xs, ys, zs = spher_2_cart(rad, theta, phi)

        ax.plot(xs, ys, zs, "b", alpha=0.5)
        


#trick to force equal unit scale
set_axes_equal(ax)
ax.view_init(elev=20., azim=-75)

plt.show()

#!/usr/bin/python

from __future__ import division

import sys
import os

import h5py


import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.patches import Ellipse
from matplotlib.colors import LogNorm 

print "Matplotlib version", matplotlib.__version__

cmap = plt.get_cmap('plasma_r')

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.ndimage as ndimage
import scipy.misc as misc
from scipy.spatial import ConvexHull, Delaunay



plt.rc('font', family='serif')
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)


#fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
fig = plt.figure(figsize=(7.48, 1.6))  #two column figure




#if len(sys.argv) != 3:
#    print sys.argv[0], ": INFILE OUTFILE"
#    exit(0)
#
#filename = sys.argv[1]
#outfilename = sys.argv[2]


#filename = 'extreme_400pts_1.4m_15km_600hz_15inc.h5'
#outfilename = 'extreme_400pts_1.4m_15km_600hz_15inc.pdf'

filename = 'casual_400pts_1.6m_12km_400hz_15inc.h5'
outfilename = 'casual_400pts_1.6m_12km_400hz_15inc.pdf'



path = os.path.abspath(os.path.dirname(sys.argv[0]))
print "Path is", path
print "Reading file", filename

f = h5py.File(filename, 'r')

# load data matrices and header variables

angular_velocity = f["header"].attrs["angular_velocity"]

redshift_matrix   = f["image/redshift_matrix"][()]
time_delay_matrix = f["image/time_delay_matrix"][()]

hit_matrix        = f["image/hit_matrix"][()]
hit_indices = hit_matrix >= 1

x_matrix          = f["image/x_matrix"][()]
y_matrix          = f["image/y_matrix"][()]
theta_matrix      = f["image/theta_matrix"][()]
phi_matrix        = f["image/phi_matrix"][()]
r_matrix          = f["image/r_matrix"][()]

angle_matrix_rest = f["image/angle_matrix_rest"][()]
angle_matrix_obs  = f["image/angle_matrix_observer"][()]

# fix nans
def fixnans(mat):
    mat[np.logical_not(np.isfinite(mat))] = 0.0
    return mat

if 1:
    fixnans(redshift_matrix)
    fixnans(time_delay_matrix)
    fixnans(theta_matrix)
    fixnans(phi_matrix)



# extent of the image plane
min_x = np.amin(x_matrix)
max_x = np.amax(x_matrix)
min_y = np.amin(y_matrix)
max_y = np.amax(y_matrix)

print 'x min/max {} {}'.format(min_x, max_x)
print 'y min/max {} {}'.format(min_y, max_y)

# in pixels
x_pixels = x_matrix.shape[0]
y_pixels = x_matrix.shape[1]

# maximum carters constant and H values 
C1_matrix = np.zeros_like(x_matrix)
C2_matrix = np.zeros_like(x_matrix)
H_matrix = np.zeros_like(x_matrix)

#
# other useful stuff

# luminosity distance
lumdist = f["header"].attrs["luminosity_distance"]

# stefan-boltzmann constant in planck units
sigma_sb = 35392.0;

# temperature conversion factor from keV to Kelvin
kelvin_per_keV = 1.16045e7;

# juris NS data
#juridata       = np.loadtxt(path + "/juridata/nu1Hz_blackbody_rho30deg.dat")
#juridata2      = np.loadtxt(path + "/juridata/nu400Hz_blackbody_rho30deg.dat")
#juridata_1deg  = np.loadtxt(path + "/juridata/nu1Hz_blackbody_rho1deg.dat")
#juridata_1deg2 = np.loadtxt(path + "/juridata/nu400Hz_blackbody_rho1deg.dat")

# imshow with some defaults
def imgplane_imshow(obj, data, interpolation="none", cmap=cmap, **kwargs):
    return obj.imshow(np.swapaxes(data, 0, 1), interpolation=interpolation, cmap=cmap, origin="lower", 
            extent=(min_x, max_x, min_y, max_y), **kwargs)


# get carter matrix data
num_geodesics = f["geodesics"].attrs["num_geodesics"]
print "geodesics in file:", num_geodesics
geodata = {}
ij_matrix = f["image/ij_matrix"]
max_redshift = np.nanmax(redshift_matrix)
print "max redshift", max_redshift

for geo_index in xrange(num_geodesics):
#for geo_index in xrange(100):
    dataset = "geodesics/geodesic_{0}/coordinates".format(geo_index)
    print "reading in geodesic data", geo_index, "from dataset", dataset

    dataset = "geodesics/geodesic_{0}/hamiltonian".format(geo_index)
    Hvalues_in = f[dataset][()]
    Hvalues = np.abs(Hvalues_in-Hvalues_in[-1])

    Cvalues_in1 = f["geodesics/geodesic_{0}/carters_constant_1".format(geo_index)][()]
    #Cvalues1 = np.abs(Cvalues_in1-Cvalues_in1[-1])
    Cvalues1 = np.abs(Cvalues_in1-Cvalues_in1[-1])/Cvalues_in1[-1]

    Cvalues_in2 = f["geodesics/geodesic_{0}/carters_constant_2".format(geo_index)][()]
    #Cvalues2 = np.abs(Cvalues_in2-Cvalues_in2[-1])
    Cvalues2 = np.abs(Cvalues_in2-Cvalues_in2[-1])/Cvalues_in2[-1]

    #Cratio = np.abs(Cvalues2/Cvalues1)

    geo_ij = ij_matrix[geo_index, :]
    H_matrix[geo_ij[0], geo_ij[1]] = np.amax(Hvalues)
    C1_matrix[geo_ij[0], geo_ij[1]] = np.amax(Cvalues1)
    C2_matrix[geo_ij[0], geo_ij[1]] = np.amax(Cvalues2)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    #return r'${} \times 10^{{{}}}$'.format(a, b)
    return r'$10^{{{}}}$'.format(b)

def plot_mat(ax, mat, title='', fmt=None, vmin=None, vmax=None, extent=None):
    if fmt is not None:
        formatter = matplotlib.ticker.FuncFormatter(fmt)
    else:
        formatter = None
    #im = ax.imshow(mat.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    im = ax.imshow(mat.T, origin='lower', cmap=cmap, 
                    norm=LogNorm(vmin=vmin, vmax=vmax,),
                    extent=extent)


    #divider = make_axes_locatable(ax)
    #cax = divider.new_horizontal(size="5%", pad=0.7, pack_start=True)
    #fig.add_axes(cax)

    fig.colorbar(im, ax=ax, shrink=0.6, format=formatter)

    ax.set_title(title)
    #ax.minorticks_on()
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    #ax.set_xlim(-10.0, 10.0)
    #ax.set_ylim(-10.0, 10.0)
    ax.set_xlim(-8.0, 8.0)
    ax.set_ylim(-8.0, 8.0)


##################################################

#Construct output xy image plane from img object
##################################################
x_span = 12.5
y_span = 12.5

x_bins = 500
y_bins = 500

xs = np.linspace(-x_span, x_span, x_bins)
ys = np.linspace(-y_span, y_span, y_bins)

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

#read redshift array
#fname = "../out/reds_f600pbbr15m1.4i15.csv"
fname = "../out/reds_f400pbbr12m1.6i15.csv"
data = np.genfromtxt(fname, delimiter=',')
redshift = np.reshape(data, (x_bins, y_bins) )
redshift = clean_image(redshift)


##################################################
#fname2 = '../out/reds_f600_bb_r15_m1.4_i15.csv'
fname2 = '../out/reds_f400_bb_r12_m1.6_i15.csv'
data2 = np.genfromtxt(fname2, delimiter=',')
redshift2 = np.reshape(data2, (x_bins, y_bins) )
redshift2 = clean_image(redshift2)


# other settings for imshow
extent=( xs[0], xs[-1], ys[0], xs[-1] )
interpolation = 'nearest'

# relative error
relerr = np.zeros(( x_bins, y_bins))
for i, x in enumerate(xs):
    for j, y in enumerate(ys):

        val1 = redshift[i,j]
        val2 = redshift2[i,j]

        errval = 0.0
        if not(val2 == 0.0):
            errval = np.abs( (val2 - val1)/val2 )
            #errval = np.log10( np.abs((val2 - val1)/val2) )
        relerr[i,j] = errval
relerr = np.ma.masked_where(relerr == 0, relerr) 
relerr = relerr.T


##################################################





gs = plt.GridSpec(1, 4)


ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[0,2])
ax4 = plt.subplot(gs[0,3])

extentC = (-8.0, 8.0, -8.0, 8.0)
#extentC = (-10.0, 10.0, -10.0, 10.0)

plot_mat(ax1, H_matrix,  title='Abs. err. in $H$', fmt=fmt,vmin=1.0e-13, vmax=1.0e-10, extent=extentC)
plot_mat(ax2, C1_matrix, title='Rel. err. in $C_1$', fmt=fmt, vmin=1.0e-3, vmax=1.0e1, extent=extentC)
plot_mat(ax3, C2_matrix, title='Rel. err. in $C_2$', fmt=fmt, vmin=1.0e-3, vmax=1.0e1, extent=extentC)
plot_mat(ax4, relerr, title='Rel. err. in $z$', fmt=fmt, vmin=1.0e-4, vmax=1.0e-2, extent=extent)


#plot_mat(ax1, H_matrix, fmt=fmt,vmin=1.0e-13, vmax=1.0e-10, extent=extentC)
#plot_mat(ax2, C1_matrix, fmt=fmt, vmin=1.0e-3, vmax=1.0e1, extent=extentC)
#plot_mat(ax3, C2_matrix, fmt=fmt, vmin=1.0e-3, vmax=1.0e1, extent=extentC)
#plot_mat(ax4, relerr, fmt=fmt, vmin=1.0e-4, vmax=1.0e-2, extent=extent)


fig.text(0.06, 0.5, '$R_{\\mathrm{e}} = 12$km \n $M=1.6 M_{\\odot}$ \n $\\nu=400$ Hz \n $i = 15^{\\circ}$', 
        ha='center', va='center', rotation='vertical' )#, size=28)
#fig.text(0.06, 0.5, '$R_{\\mathrm{e}} = 15$km \n $M=1.4 M_{\\odot}$ \n $\\nu=600$ Hz \n $i = 15^{\\circ}$', 
#        ha='center', va='center', rotation='vertical' )#, size=28)


plt.subplots_adjust(left=0.12, bottom=-0.13, right=0.98, top=1, wspace=0.11, hspace=0)

#fig.tight_layout()
#plt.tight_layout()
#plt.subplots_adjust(wspace=0.3)
#plt.savefig(outfilename )#, bbox_inches='tight')
plt.savefig(outfilename, bbox_inches='tight')
#plt.show()



from __future__ import division

import sys
import os


import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

print "Matplotlib version", matplotlib.__version__

from matplotlib import cm
import palettable as pal


#cmap = pal.cmocean.sequential.Matter_8.mpl_colormap #best so far
#cmap = pal.wesanderson.Zissou_5.mpl_colormap
cmap = pal.colorbrewer.qualitative.Set1_6.mpl_colormap
#cmap = plt.get_cmap('plasma_r')
#cmap = cm.get_cmap('inferno_r')

import numpy as np
from scipy.signal import savgol_filter



plt.rc('font', family='serif')
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)


#fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
fig = plt.figure(figsize=(7.48, 4.0))  #two column figure

outfilename = 'sweep_grid.pdf'


#Read JN files
def read_csv(fname):
    da = np.genfromtxt(fname, delimiter=",")

    des = np.diff(da[:,0])[2]
    norm = np.sum(des*da[:,1])

    return da[:,0],da[:,1] #/norm


def kernel_smooth(xx, yy):
    kernel = [ 0.1174,  0.1975,  0.2349,  0.1975,  0.1174 ]
    for i in range(3, len(xx)-3):
        for j,indx in enumerate( [-2,-1,0,1,2] ):
            yy[i] += kernel[j] * yy[j + indx]
    return xx, yy

def smooth(xx, yy):
    yy = savgol_filter(yy, 5, 1)
    np.clip(yy, 0.0, 1000.0, out=yy)
    yy[0] = 0.0
    yy[-1] = 0.0
    return xx, yy



def fwhm(xx, yy):

    lh = 0.0
    rh = 0.0
    lf = 0.0
    rf = 0.0

    yimax = np.argmax(yy)
    for i, x in enumerate(xx):
        if i < yimax:
            if yy[i] < 0.01:
                lf = x
            if yy[i] < 0.5:
                lh = x
        else:
            if yy[i] >= 0.01:
                rf = x
            if yy[i] >= 0.5:
                rh = x

    return lh, rh, lf, rf


solar_mass_per_km = 0.6772

freq = 600

xmin = 0.9
xmax = 1.05

xmin = 0.76
xmax = 1.28

#gs = plt.GridSpec(3, 7)
gs = plt.GridSpec(3, 4)

for i, mass in enumerate([1.8, 1.5, 1.1]):

    #for j, rad in enumerate([10,11,12,13,14,15,16]):
    for j, rad in enumerate([10,12,14,16]):

        ax = plt.subplot(gs[i, j])
        ax.minorticks_on()
        ax.set_ylim(0.0, 1.1)
        ax.set_xlim(xmin, xmax)

        if i == 0:
            string = r'$R$=' + str(rad) + ' km'
            ax.set_title(string, ha='center', va='center', size=12)
        if i == 2:
            ax.set_xlabel(r'Relative energy $E/\sqrt{ 1-u }$')

        if not( j == 0 ):
            ax.set_yticklabels( [] )

        if (i == 1) and (j == 0):
            ax.set_ylabel(r'Peak normalized flux')

        if j == 3:
            axl = ax.twinx()
            axl.set_yticklabels( [] )
            axl.set_ylabel(r'$M=$'+str(mass), size=12, rotation=270, labelpad=10)



        radius = rad * solar_mass_per_km / mass
        compactness = np.sqrt(1 - 2/radius)



        if not( i == 2 ):
            ax.set_xticklabels( [] )
        #else:
        #    plt.setp( ax.xaxis.get_majorticklabels(), rotation=60 )


        ################################################## 
        #now plot

        #col = 'black'
        #for k, incl in enumerate([90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 5, 1]):
        #for k, incl in enumerate([40, 30, 20, 15, 10, 5, 1]):
        #for k, incl in enumerate([1, 5, 10, 15, 20, 40, 60, 90]):
        for k, incl in enumerate([90, 60, 40, 20, 10, 5]):
            fname = 'sweep/lineprofile_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.csv'.format(
            np.int(freq), np.int(rad), mass, np.int(incl))

            
            if not( os.path.isfile( fname ) ):
                continue

            col = cmap( (5-k) / float(5) )
            if k == 0:
                col = 'black'

            xx, yy = read_csv(fname)

            xx, yy = smooth(xx, yy)

            yy /= np.max(yy)
            xx /= compactness

            lh, rh, lf, rf = fwhm(xx, yy)

            ax.plot( xx, yy, linestyle='solid', color=col )

            #ax.plot([lh, lh], [0.0, 1.0], color=col, linestyle='solid')
            #ax.plot([rh, rh], [0.0, 1.0], color=col, linestyle='solid')

            #ax.plot([lf, lf], [0.0, 1.0], color=col, linestyle='solid')
            #ax.plot([rf, rf], [0.0, 1.0], color=col, linestyle='solid')







plt.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
plt.savefig(outfilename, bbox_inches='tight')

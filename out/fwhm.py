from __future__ import division

import sys
import os


import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

print "Matplotlib version", matplotlib.__version__

from matplotlib import cm
import palettable as pal


from scipy.signal import savgol_filter
import numpy as np


plt.rc('font', family='serif')
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
#plt.rcParams['image.cmap'] = 'inferno_r'
plt.rcParams['image.cmap'] = 'viridis_r'


#fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
fig = plt.figure(figsize=(7.48, 3.0))  #two column figure


#mode = 'abs'
mode = 'rel'
#mode = 'full'

if mode == 'abs':
    label = 'FWHM'
elif mode == 'rel':
    label= 'FWHM / FWTM'
elif mode == 'full':
    label = 'FWTM'
fig.text(0.98, 0.5, label, size = 13, ha='center', va='center', rotation=270)



outfilename = 'fwhm_'+mode+'.pdf'


#Read JN files
def read_csv(fname):
    da = np.genfromtxt(fname, delimiter=",")

    des = np.diff(da[:,0])[2]
    norm = np.sum(des*da[:,1])

    return da[:,0],da[:,1] #/norm

def smooth(xx, yy):
    yy = savgol_filter(yy, 5, 2)
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
            if yy[i] < 0.1:
                lf = x
            if yy[i] < 0.5:
                lh = x
        else:
            if yy[i] >= 0.1:
                rf = x
            if yy[i] >= 0.5:
                rh = x

    return lh, rh, lf, rf




solar_mass_per_km = 0.6772

freq = 600

xmin = 10.0
xmax = 16.0
ymin = 0.0
ymax = 90.0


gs = plt.GridSpec(1, 3)

mass_grid = [1.1, 1.5, 1.8]
rad_grid = [10,11,12,13,14,15,16]
incl_grid = [90,80,70,60,50,40,30,20,15,10,9,8,7,6,5,4,3,2,1,0.5]

for i, mass in enumerate(mass_grid):

    xg1 = np.zeros(len(rad_grid) * len(incl_grid))
    yg1 = np.zeros(len(rad_grid) * len(incl_grid))
    zg1 = np.zeros(len(rad_grid) * len(incl_grid))

    grid = np.zeros( (len(rad_grid), len(incl_grid)) )

    ax = plt.subplot(gs[0, i])
    ax.minorticks_on()

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    if mode == 'full':
        string = r'$M=$'+str(mass)+'M$_{\odot}$'
        ax.set_title(string, ha='center', va='center', size=12)

    if mode == 'rel':
        ax.set_xlabel( r'Radius $R_{\mathrm{e}}$ (km)', size=8)

    if not( i == 0 ):
        ax.set_yticklabels( [] )
    else:
        ax.set_ylabel(r'Inclination $i$ (deg)', size=8)

    q = 0
    for j, rad in enumerate(rad_grid):


        radius = rad * solar_mass_per_km / mass
        compactness = np.sqrt(1 - 2/radius)

        for k, incl in enumerate(incl_grid):

            fname = 'sweep/lineprofile_f{:03d}_bb_r{:02d}_m{:03.1f}_i{:02d}.csv'.format(
            np.int(freq), np.int(rad), mass, np.int(incl))

            
            if not( os.path.isfile( fname ) ):
                continue

            xx, yy = read_csv(fname)
            #xx, yy = smooth(xx, yy)
            yy /= np.max(yy)
            #xx /= compactness

            lh, rh, lf, rf = fwhm(xx, yy)
            
            fullw = np.abs(rf-lf)
            halfw = np.abs(rh-lh)

            if mode == 'abs':
                grid[j, k] = halfw
                zg1[q] = halfw
            elif mode == 'rel':
                grid[j, k] = halfw/fullw
                zg1[q] = halfw/fullw
            elif mode == 'full':
                grid[j, k] = fullw
                zg1[q] = fullw

            #grid[j, k] = fullw
            #grid[j, k] = halfw/fullw
            #grid[j, k] = fullw

            xg1[q] = rad
            yg1[q] = incl
            q += 1


    ##################################################
    if mode == 'abs':
        vmin = 0.0
        vmax = 0.5
    elif mode == 'rel':
        vmin = 0.0
        vmax = 1.0
    elif mode == 'full':
        vmin = 0.0
        vmax = 0.42

    xg, yg = np.meshgrid(rad_grid, incl_grid)


    interpolation = 'bicubic'

    levels = np.linspace(vmin, vmax, 100)

    xg1 = xg1[0:q]
    yg1 = yg1[0:q]
    zg1 = zg1[0:q]

    if False:
        #ax.tripcolor(xg,yg,zg)
        CS = ax.tricontourf(xg1,yg1,zg1, levels)
        #ax.plot(xg1,yg1, 'k.')


    if True:
        CS = ax.contourf(xg, yg, grid.T,
                   levels,
                   #colors = 'w',
                   linestyle = 'dashed',
                   interpolation=interpolation,
                   vmin = vmin,
                   vmax = vmax,
                   )


    # This is the fix for the white lines between contour levels
    for c in CS.collections:
        c.set_edgecolor("face")

    #cb = plt.colorbar(cnt, ticks=np.linspace(vmin, vmax, 5) )
    #levels = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    #ax.contour(xg, yg, grid.T,
    #           levels,
    #           colors = 'w',
    #           linestyle = 'dashed',
    #           interpolation='bicubic',
    #           vmin = vmin,
    #           vmax = vmax,
    #           )

    if mode == 'abs':
        levels = [0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4]
        def fmt(x):
            return r'$\Delta E/E$ = {:.2f}'.format(x)

    elif mode == 'rel':
        #levels = [0.1, 0.2, 0.4, 0.6, 0.8, 0.85, 0.9, 0.95]
        levels = [0.1, 0.5, 0.8, 0.85, 0.9, 0.95]
        def fmt(x):
            a = 100*x
            return '{:2d}%'.format(np.int(a))

    elif mode == 'full':
        levels = [0.01, 0.05, 0.1, 0.2, 0.3]
        def fmt(x):
            return r'$\Delta E/E$ = {:.2f}'.format(x)

    if False:
        CS = ax.tricontour(xg1,yg1,zg1, 
                levels,
                colors = 'k',
                linestyle='solid',
                )


    if True:
        CS = ax.contour(xg, yg, grid.T,
                   levels,
                   colors = 'k',
                   linestyle = 'solid',
                   interpolation=interpolation,
                   vmin = vmin,
                   vmax = vmax,
                   )

    plt.clabel(CS, inline=1, fontsize=6, inline_spacing=1, fmt=fmt)

plt.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.85, wspace=0.1, hspace=0.1)
#plt.savefig(outfilename )#, bbox_inches='tight')
plt.savefig(outfilename, bbox_inches='tight')

import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

import scipy.ndimage as ndimage

#read JP and TH files
#def read_JP_files(fname):
#    da = np.genfromtxt(fname, delimiter="    ", comments='#')
#    return da[:,0], da[:,1], da[:,2], da[:,3],da[:,4],da[:,5]

#Read JN files
def read_JN_files(fname):
    da = np.genfromtxt(fname, delimiter=",")
    return da[:,0],da[:,1],da[:,2],da[:,3],da[:,4],da[:,5],da[:,6],da[:,7],da[:,8]


    
## Plot
fig = figure(figsize=(9,4), dpi=80)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 100)
gs.update(wspace = 0.34)
#gs.update(hspace = 0.4)


lsize = 7.0

#phase limits
xmin = 0.0
xmax = 90.0

ymin = 0.0
ymax = 90.0



#figure shape parameters
panelh = 70
skiph = 30

mfiglim = 0

#path to files
path_files = "../../out_skymaps/"

#labels size
tsize = 8.0

#general parameters
nu = 'f600'
bprof = 'pbb'
rad = 'r15'
mass = 'm1.6'
rho = 'x10'

incls = ['i5','i10','i20','i30','i40','i50','i60','i70','i80','i90']
incls_g = [5,10,20,30,40,50,60,70,80,90]

colat_g = [10,30,50,70,90]


#pre-read one file to get initial values
colat = 'd10'
incl = incls[0]
fname = path_files + nu+bprof+rad+mass+colat+incl+rho
phase_g, N2kev, N6kev, N12kev, Nbol, Fbol, F2kev, F6kev, F12kev = read_JN_files(fname+'.csv') 
Nt = len(phase_g)


phase_t = np.linspace(0.0, 1.0,  200)
incls_t = np.linspace(0.0, 90.0, 100)

maxflux = 0.0
fig.text(0.3, 0.92, 'North pole spot',  ha='center', va='center', size=10)
fig.text(0.7, 0.92, 'North + south pole spots',  ha='center', va='center', size=10)


#empty matrix
pfracs = np.zeros((len(colat_g), len(incls_g)))


for k in range(2):

    #frame for the main pulse profile fig
    #ax1 = subplot(gs[mfiglim:mfiglim+panelh, k])
    if k == 0:
        #ax1 = subplot(gs[mfiglim:mfiglim+panelh, 0:46])
        ax1 = subplot(gs[0, 0:46])
    else:
        #ax1 = subplot(gs[mfiglim:mfiglim+panelh, 49:95])
        ax1 = subplot(gs[0, 49:95])

    ax1.minorticks_on()
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    
    ax1.set_xlabel('Inclination $i$', size=lsize)

    if k == 0:
        ax1.set_ylabel('Spot colatitude $\\theta_s$', size=lsize)
    elif k == 1:
        ax1.set_yticklabels([])

        
    for j in range(5):
                
        if j == 0:
            colat = '10'
        elif j == 1:
            colat = '30'
        elif j == 2:
            colat = '50'
        elif j == 3:
            colat = '70'
        elif j == 4:
            colat = '90'
        
            
        #skymap = np.zeros((Nt, len(incls)))
        #skymap = np.zeros((len(incls), Nt))
        
        for q in range(len(incls)):
            incl = incls[q]
            #incl = incls[0]
            
            fname = path_files + nu+bprof+rad+mass+'d'+colat+incl+rho
            phase, N2kev, N6kev, N12kev, Nbol, Fbol, F2kev, F6kev, F12kev = read_JN_files(fname+'.csv') 

            #add second spot
            if k == 1:
                phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2, F2kev2, F6kev2, F12kev2 = read_JN_files(fname+'_2nd.csv') 
                N2kev  += N2kev2 
                N6kev  += N6kev2
                N12kev += N12kev2
                Nbol   += Nbol2
                Fbol   += Fbol2
                F2kev  += F2kev2
                F6kev  += F6kev2
                F12kev += F12kev2

            #build flux matrix
            flux = Fbol
            #flux = Fbol / flux.max()

            #skymap[:,q] = flux
            #skymap[q,:] = flux

            pfracs[j,q] = (flux.max() - flux.min()) / (flux.max() + flux.min())
            
    #print skymap.max()
    #print shape(skymap)
    #print skymap
    #skymap_interp = griddata((phase_g, incls_g), skymap, (phase_t, incls_t), method='cubic')
    #skymap_interp = griddata((phase_g, incls_g), skymap, np.meshgrid(phase_t, incls_t), method='cubic')
    #print skymap_interp

    hdata = pfracs
    
    xr0 = incls_g[0]
    xr1 = incls_g[-1]
    yr0 = colat_g[0]
    yr1 = colat_g[-1]


    #print xr0, xr1, yr0, yr1
    extent = [xr0, xr1, yr0, yr1]

    
    hdata_smooth = ndimage.gaussian_filter(hdata, sigma=1.0, order=0)
    #hdata_masked = np.ma.masked_where(hdata <= 0.001, hdata)
    #im = ax1.imshow(hdata_masked.T,
    im = ax1.imshow(hdata.T,
                    #interpolation='nearest',
                    interpolation='gaussian',
                    origin='lower',
                    extent=extent,
                    #cmap='Reds',
                    #cmap='jet',
                    #cmap='YlGnBu',
                    cmap='plasma_r',
                    vmin=0.0,
                    #vmax=0.4,
                    vmax=1.0,
                    aspect='auto')

    levels = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    #levels = [0.05, 0.1, 0.2, 0.25, 0.3, 0.35]
    cs1 = ax1.contour(hdata_smooth.T,
    #cs1 = ax1.contour(hdata.T,
                      levels,
                      colors = 'r',
                      origin='lower',
                      extent=extent)
    zc = cs1.collections[0]
    setp(zc, linewidth=1)

    print hdata
        
    if k == 1:
        #mfiglim:mfiglim+panelh, 0:40])
        #cbaxes = fig.add_axes([0.90, (mfiglim+panelh)/500, 0.05, panelh/500.0])
        cbaxes = subplot(gs[0, 95:97])
        cb = colorbar(im,
                      #label='Probability density',
                      cax=cbaxes)
        cb.set_label('Pulse fraction',size=lsize)


    #fig.text(0.5, 0.91-j*0.16, '$\\theta_{\mathrm{s}}$ = '+colat,  ha='center', va='center', size=tsize)
    mfiglim += panelh+skiph

    


savefig('fig8.pdf', bbox_inches='tight')

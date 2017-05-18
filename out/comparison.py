import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

#from sg_filter import savitzky_golay
from scipy.signal import savgol_filter
import scipy.fftpack

def read_JP_files(fname):
    da = np.genfromtxt(fname, delimiter="   ")
    return da[:,0], da[:,1], da[:,2], da[:,3],da[:,4],da[:,5]

def read_JN_files(fname):
    da = np.genfromtxt(fname, delimiter=",")
    return da[:,0],da[:,1],da[:,2],da[:,3],da[:,4],da[:,5]

def read_PP_files(fname):
    da = np.genfromtxt(fname, delimiter=" ")
    return da[:,0],da[:,1]


    
## Plot
fig = figure(figsize=(9,10), dpi=80)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(400, 4)
gs.update(wspace = 0.34)
#gs.update(hspace = 0.4)


lsize = 7.0

xmin = -0.04
xmax = 1.04


#eymin = -1.0
#eymax = 1.0
eymin = -0.5
eymax = 0.5


panelh = 45
epanelh = 25
skiph = 30


mfiglim = 0

path_JP = ""

#labels
tsize = 10.0


#nu = '1'
nu = '400'
fig.text(0.5, 0.92, '$\\nu = 1  $ Hz  blackbody  $\\rho = 30^{\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.72, '$\\nu = 400$ Hz  blackbody  $\\rho = 30^{\circ}$',  ha='center', va='center', size=tsize)

#fig.text(0.5, 0.52, '$\\nu = '+nu+'$ Hz  blackbody  $\\rho = 30^{\circ}$',  ha='center', va='center', size=tsize)
#fig.text(0.5, 0.32, '$\\nu = '+nu+'$ Hz  Hopf  $\\rho = 30^{\circ}$',  ha='center', va='center', size=tsize)

#fig.text(0.5, 0.12, 'Phase',ha='center', va='center', size=lsize)


for j in range(2):

    if j == 0:
        fname  = path_JP + 'nu1Hz_blackbody_rho30deg.dat'
        fname2 = path_JP + 'f001_bb_r12_m1.6_d50_i60_x30.csv'
    if j == 1:
        fname  = path_JP + 'nu400Hz_blackbody_rho30deg.dat'
        fname2 = path_JP + 'f400_bb_r12_m1.6_d50_i60_x30.csv'
    




    #read JP data
    phase, N2kev, N6kev, N12kev, Nbol, Fbol = read_JP_files(fname)

    #read JN data
    phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2 = read_JN_files(fname2) 
    phase2s = phase2

    
    for i in range(4):


         #frame for the main pulse profile fig
         ax1 = subplot(gs[mfiglim:mfiglim+panelh, i])
         ax1.minorticks_on()
         #ax1.set_xticklabels([])
         ax1.set_xlim(xmin, xmax)

         #ax1.yaxis.major.formatter.set_powerlimits((0,0))

         formatter = ScalarFormatter(useMathText=True)
         formatter.set_scientific(True)
         formatter.set_powerlimits((0,0))
         ax1.yaxis.set_major_formatter(formatter)


         fscale = 1.0
         
         if i == 0:
             ax1.set_ylabel('$N$ (2 keV)\n[ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$]',size=lsize)
             flux = N2kev
             flux2 = N2kev2/fscale
         elif i == 1:
             ax1.set_ylabel('$N$ (6 keV)',size=lsize)
             flux = N6kev
             flux2 = N6kev2/fscale
         elif i == 2:
             ax1.set_ylabel('$N$ (12 keV)',size=lsize)
             flux = N12kev
             flux2 = N12kev2/fscale
         elif i == 3:
             ax1.set_ylabel('Bolometric [ph cm$^{-2}$ s$^{-1}$]',size=lsize)
             flux = Nbol
             flux2 = Nbol2/fscale
             #flux  = Fbol
             #flux2 = Fbol2
             


         #Savitzky-Golay low-pass filtter
         flux2[-1] = flux2[0]
         #flux3 = savgol_filter(flux2, 15, 5, mode='wrap')
         flux3 = savgol_filter(flux2, 15, 5, mode='wrap')
         flux3[0:9] = flux2[0:9]
         flux3[-9:-1] = flux2[-9:-1]
         flux3[-1] = flux2[0]

         #JP data
         ax1.plot(phase, flux, 'k-')

         #JN data
         ax1.plot(phase2, flux2, 'r--')
             

         ax1.set_yticks(ax1.get_yticks()[1:-1])
         
         #frame for the error panel
         ax2 = subplot(gs[(mfiglim+panelh):(mfiglim+panelh+epanelh), i])
         ax2.minorticks_on()
         ax2.set_xlim(xmin, xmax)
         #ax2.set_ylim(eymin, eymax)
         #ax2.set_yticks(ax2.get_yticks()[1:-1])

         if i == 0:
             ax2.set_ylabel('$\Delta$ %',size=lsize)
             
         if j == 3:
            ax2.set_xlabel('Phase', size=lsize)
            
         ax2.plot([xmin, xmax], [0.0, 0.0], 'r--', linewidth=0.3)


         #interpolate error from JN
         #fluxi2 = griddata(phase2, flux2, (phase), method='cubic')
         #fluxi2 = griddata(phase2, flux2, (phase), method='linear')
         #err = (flux/fluxi2 - 1)*100
         #ax2.plot(phase, err, 'k-', linewidth = 0.4)
         
         #interpolate error from JP
         #fluxi = griddata(phase, flux, (phase2), method='linear')
         fluxi = griddata(phase, flux, (phase2), method='cubic')
         err = (fluxi/flux2 - 1)*100
         #ax2.plot(phase2, err, 'k-', linewidth = 0.4)

         err = (fluxi/flux3 - 1)*100
         ax2.plot(phase2, err, 'k-', linewidth = 0.4)


    mfiglim += panelh+epanelh+skiph

    


savefig('pulse.pdf', bbox_inches='tight')

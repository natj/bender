import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

def read_JP_files(fname):
    da = np.genfromtxt(fname, delimiter="    ", comments='#')
    return da[:,0], da[:,1], da[:,2], da[:,3],da[:,4],da[:,5]

def read_JN_files(fname):
    da = np.genfromtxt(fname, delimiter=",")
    return da[:,0],da[:,1],da[:,2],da[:,3],da[:,4],da[:,5],da[:,6],da[:,7],da[:,8]


    
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


eymin = -5.0
eymax = 5.0


panelh = 45
epanelh = 25
skiph = 30


mfiglim = 0

path_JP = "../../out2/f700/r12nn/"

#labels
tsize = 10.0


nu = '700'
#nu = '400'


fig.text(0.5, 0.92, '$\\theta_s = 18^{\\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.72, '$\\theta_s = 45^{\\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.52, '$\\theta_s = 90^{\\circ}$',  ha='center', va='center', size=tsize)

#fig.text(0.5, 0.12, 'Phase',ha='center', va='center', size=lsize)


for j in range(3):

    if j == 0:
        fname = path_JP + 'r12x10d18i45.txt'
        fname2 = path_JP + 'f'+nu+'pbbr12m1.4d18i45x10.csv'
    if j == 1:
        fname = path_JP + 'r12x10d45i45.txt'
        fname2 = path_JP + 'f'+nu+'pbbr12m1.4d45i45x10.csv'
    if j == 2:
        fname = path_JP + 'r12x10d90i45.txt'
        fname2 = path_JP + 'f'+nu+'pbbr12m1.4d90i45x10.csv'
    #if j == 3:
    #    fname = path_JP + 'nu'+nu+'Hz_hopf_rho30deg.dat'
    #    fname2 = path_JP + 'f'+nu+'phopfr12m1.6d50i60x30.csv'
    

    #read JP data
    phase, N2kev, N6kev, N12kev, Nbol, Fbol = read_JP_files(fname)

    #read JN data
    phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2, F2kev2, F6kev2, F12kev2 = read_JN_files(fname2) 

    
    
    for i in range(4):


         #frame for the main pulse profile fig
         ax1 = subplot(gs[mfiglim:mfiglim+panelh, i])
         ax1.minorticks_on()
         ax1.set_xticklabels([])
         ax1.set_xlim(xmin, xmax)

         #ax1.yaxis.major.formatter.set_powerlimits((0,0))

         formatter = ScalarFormatter(useMathText=True)
         formatter.set_scientific(True)
         formatter.set_powerlimits((0,0))
         ax1.yaxis.set_major_formatter(formatter)
         
         #ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

         #xmft = ScalarFormatter()
         #xmft.set_powerlimits((-2,2))
         #ax1.xaxis.set_major_formatter(xmft)
         
         if i == 0:
             ax1.set_ylabel('$N$ (2 keV)\n[ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$]',size=lsize)
             flux = N2kev
             flux2 = N2kev2
         elif i == 1:
             ax1.set_ylabel('$N$ (6 keV)',size=lsize)
             flux = N6kev
             flux2 = N6kev2
         elif i == 2:
             ax1.set_ylabel('$N$ (12 keV)',size=lsize)
             flux = N12kev
             flux2 = N12kev2
         elif i == 3:
             ax1.set_ylabel('Bolometric [ph cm$^{-2}$ s$^{-1}$]',size=lsize)
             flux = Nbol
             flux2 = Nbol2

             
         #JP data
         ax1.plot(phase, flux, 'k-')

         #flux2 = flux2 * 0.99
         #phase2 = phase2 - 0.0018
         phase2 = phase2 - 0.0027
         #phase = phase - 0.01
         
         #JN data
         ax1.plot(phase2, flux2, 'r:')
         #ax1.plot(phase2, flux2, 'r-', linewidth=0.3)
         
         #frame for the error panel
         ax2 = subplot(gs[(mfiglim+panelh):(mfiglim+panelh+epanelh), i])
         ax2.minorticks_on()
         ax2.set_xlim(xmin, xmax)
         ax2.set_ylim(eymin, eymax)

         if i == 0:
             ax2.set_ylabel('$\Delta$ %',size=lsize)

             
         #if j != 3:
         #   ax2.set_xticklabels([])

         if j == 2:
            ax2.set_xlabel('Phase', size=lsize)
            
         ax2.plot([xmin, xmax], [0.0, 0.0], 'r--', linewidth=0.3)


         #interpolate error
         #fluxi = interp1d(phase, flux, kind='linear')
         fluxi2 = griddata(phase2, flux2, (phase), method='cubic')
         
         #fluxi = interp1d(phase, flux, kind='cubic')
         err = (fluxi2/flux - 1)*100
                  
         #flux2i = interp1d(phase2, flux2, kind='cubic', fill_value='extrapolate')
         #err = (flux/flux2i(phase) - 1)*100

         ax2.plot(phase, err, 'k-', linewidth = 0.4)


         for pshift in np.linspace(-0.01, 0.01, 10):
             fluxi2 = griddata(phase2+pshift, flux2, (phase), method='cubic')
             err = (fluxi2/flux - 1)*100
             #ax2.plot(phase, err, 'b-', linewidth = 0.4)

         
    mfiglim += panelh+epanelh+skiph

    


savefig('fig4.pdf', bbox_inches='tight')
#savefig('fig2b.pdf', bbox_inches='tight')

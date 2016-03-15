import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

#read JP and TH files
def read_JP_files(fname):
    da = np.genfromtxt(fname, delimiter="    ", comments='#')
    return da[:,0], da[:,1], da[:,2], da[:,3],da[:,4],da[:,5]

#Read JN files
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

#phase limits
xmin = -0.04
xmax = 1.04

#error window limits
eymin = -2.0
eymax = 2.0

#figure shape parameters
panelh = 45
epanelh = 25
skiph = 30


mfiglim = 0

#path to files
path_JP = "../../out2/f600c/"

#labels size
tsize = 10.0


nu = '600'


fig.text(0.5, 0.92, 'obl+nodel',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.72, 'obl+nodel $\\times$ $\\delta$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.52, 'obl $\\times$ $\\gamma$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.32, 'sphere $\\times$ $\\gamma$',  ha='center', va='center', size=tsize)

#fig.text(0.5, 0.12, 'Phase',ha='center', va='center', size=lsize)


for j in range(4):

    if j == 0:
        #fname = path_JP + 'rho1d50obl.txt'
        #fname = path_JP + 'rho1d50sphere.txt'
        fname = path_JP + 'rho30d50oblnodel.txt' ###
        fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_obl_nodel.csv'
    elif j == 1:
        #fname = path_JP + 'rho1d50obl.txt'
        #fname = path_JP + 'rho1d50sphere.txt'
        fname = path_JP + 'rho30d50oblnodel.txt' ###
        fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_obl_nodel_delta.csv'
    elif j == 2:
        fname = path_JP + 'rho30d50obl.txt'
        #fname = path_JP + 'rho1d50sphere.txt'
        #fname = path_JP + 'rho1d50oblnodel.txt' ###
        #fname = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_obl_gamma.csv'
        #fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_obl_gamma_polar.csv'
        #fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_obl_gamma_hires_polar.csv'
        #fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_obl_henon4.csv'
        fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_obl_test.csv'
        #fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30.csv'
    elif j == 3:
        #fname = path_JP + 'rho1d50obl.txt'
        fname = path_JP + 'rho30d50sphere.txt'
        #fname = path_JP + 'rho1d50oblnodel.txt' ###
        #fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_sphere_gamma.csv'
        #fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_sphere_gamma_polar3.csv'

        #fname = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_sphere_henon.csv'
        fname2 = path_JP + 'f'+nu+'pbbr15m1.6d50i60x30_sphere_nohenon.csv'
        


        
    #read JP data
    #if j != 2:
    phase, N2kev, N6kev, N12kev, Nbol, Fbol = read_JP_files(fname)
    #else:
    #    phase, N2kev, N6kev, N12kev, Nbol, Fbol, F2kev, F6kev, F12kev = read_JN_files(fname) 
    #read JP data
    #phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2 = read_JP_files(fname2)

    #read JN data
    #if j == 2:
    #    phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2 = read_JP_files(fname2)
    #else:
    phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2, F2kev2, F6kev2, F12kev2 = read_JN_files(fname2) 

    phasetmp = phase2
    
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
             #flux = Fbol
             #flux2 = Fbol2
             
         indxs = []
         for q in range(len(flux2)):
             if not (np.isnan(flux2[q])):
             #if not (flux2[q] == flux2[q]):
                 indxs.append(q)

         phase2 = phasetmp[indxs]
         flux2 = flux2[indxs]

         
         #JP data
         ax1.plot(phase, flux, 'k-')

         if i == 0:
             pshft = 0.0
             merr = 1.0e6
             for pshift in np.linspace(-0.1, 0.1, 100):
                 fluxi2 = griddata(phase2 + pshift, flux2, (phase), method='cubic', fill_value=0.0)
                 err = (fluxi2/flux - 1)*100

                 serr = 0.0
                 for ijk in range(len(err)):
                     if fluxi2[ijk] != 0:
                         serr += np.abs(err[ijk])
                 if serr < merr:
                     merr = serr
                     pshft = pshift
             
             print "min shift:", pshft




         
         #arbitrary phase shifts
         #flux2 = flux2 * 0.99
         #phase2 = phase2 + pshft
         if j == 0:
             phase2 = phase2 - 0.0008
         elif j == 1:
             phase2 = phase2 - 0.0008
         elif j == 2:
             phase2 = phase2 + 0.005
         elif j == 3:
             phase2 = phase2 #- 0.0005
             
         #phase = phase - 0.01
         
         #JN data
         #ax1.plot(phase2, flux2, 'r:')
         ax1.plot(phase2, flux2, 'r--')
         #ax1.plot(phase2, flux2, 'r-', linewidth=0.3)
         
         #frame for the error panel
         ax2 = subplot(gs[(mfiglim+panelh):(mfiglim+panelh+epanelh), i])
         ax2.minorticks_on()
         ax2.set_xlim(xmin, xmax)
         #ax2.set_ylim(eymin, eymax)

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
         #fluxi2 = griddata(phase2, flux2, (phase), method='linear')
         
         #fluxi = interp1d(phase, flux, kind='cubic')
         err = (flux/fluxi2 - 1)*100
                  
         #flux2i = interp1d(phase2, flux2, kind='cubic', fill_value='extrapolate')
         #err = (flux/flux2i(phase) - 1)*100

         #for q in range(len(phase)):
         #    print phase[q], err[q], fluxi2[q], flux[q]
             
         ax2.plot(phase, err, 'k-', linewidth = 0.4)


         #optional errors for range of phase shifts
         for pshift in np.linspace(-0.005, 0.005, 10):
             fluxi2 = griddata(phase2+pshift, flux2, (phase), method='cubic')
             err = (flux/fluxi2 - 1)*100
             #ax2.plot(phase, err, 'b-', linewidth = 0.4)

         
    mfiglim += panelh+epanelh+skiph

    


savefig('fig4c.pdf', bbox_inches='tight')

import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

from scipy.signal import savgol_filter



#read JP and TH files
def read_JP_files(fname):
    da = np.genfromtxt(fname, delimiter="   ", comments='#')
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
eymin = -0.5
eymax = 0.5

#figure shape parameters
panelh = 45
epanelh = 25
skiph = 30


mfiglim = 0

#path to files
#path_JP = "../../out2/f700/r12nnn/"
#path_JP = "../../out2/f700n/r12/"
path_JP = "../../out2/f700/r12/"

#labels size
tsize = 10.0


nu = '700'


fig.text(0.5, 0.92, '$\\theta_s = 18^{\\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.72, '$\\theta_s = 45^{\\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.52, '$\\theta_s = 90^{\\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.32, 'Hopf $\\theta_s = 45^{\circ}$',  ha='center', va='center', size=tsize)

#fig.text(0.5, 0.12, 'Phase',ha='center', va='center', size=lsize)


for j in range(4):

    if j == 0:
        fname = path_JP + 'iso18.txt'
        fname2 = path_JP + 'f'+nu+'pbbr12m1.4d18i45x10.csv'
    if j == 1:
        fname = path_JP + 'iso45.txt'
        fname2 = path_JP + 'f'+nu+'pbbr12m1.4d45i45x10.csv'
    if j == 2:
        fname = path_JP + 'iso90++.txt'
        fname2 = path_JP + 'f'+nu+'pbbr12m1.4d90i45x10.csv'
    if j == 3:
        #fname = path_JP + 'hopf90.txt'
        #fname2 = path_JP + 'f'+nu+'phopfr12m1.4d90i45x10_gdS.csv'
        #fname2 = path_JP + 'f'+nu+'phopfr12m1.4d90i45x10_exact.csv'
        #fname2 = path_JP + 'f'+nu+'phopfr12m1.4d90i45x10.csv'

        fname = path_JP + 'hopf45.txt'
        fname2 = path_JP + 'f'+nu+'phopfr12m1.4d45i45x10.csv'

        #fname = path_JP + 'iso90_rho3.txt'
        #fname2 = path_JP + 'f'+nu+'pbbr12m1.4d90i45x3.csv'

        #fname = path_JP + 'iso90_rho30.txt'
        #fname2 = path_JP + 'f'+nu+'pbbr12m1.4d90i45x30.csv'
        #fname3 = path_JP + 'f'+nu+'pbbr12m1.4d90i45x30_angle.csv'
    

    #read JP data
    phase, N2kev, N6kev, N12kev, Nbol, Fbol = read_JP_files(fname)

    #read JN data
    phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2, F2kev2, F6kev2, F12kev2 = read_JN_files(fname2) 

    phasetmp = phase2
    
    #read JN data reference
    #phase3, N2kev3, N6kev3, N12kev3, Nbol3, Fbol3, F2kev3, F6kev3, F12kev3 = read_JN_files(fname3) 
    #phasetmp3 = phase3

    for i in range(4):


         #frame for the main pulse profile fig
         ax1 = subplot(gs[mfiglim:mfiglim+panelh, i])
         ax1.minorticks_on()
         ax1.set_xlim(xmin, xmax)
         ax1.set_xticklabels([])
         #ax1.set_yscale('log')

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
             #fluxr = N2kev3
         elif i == 1:
             ax1.set_ylabel('$N$ (6 keV)',size=lsize)
             flux = N6kev
             flux2 = N6kev2
             #fluxr = N6kev3
         elif i == 2:
             ax1.set_ylabel('$N$ (12 keV)',size=lsize)
             flux = N12kev
             flux2 = N12kev2
             #fluxr = N12kev3
         elif i == 3:
             ax1.set_ylabel('Bolometric [ph cm$^{-2}$ s$^{-1}$]',size=lsize)
             flux = Nbol
             flux2 = Nbol2

             
         indxs = []
         for q in range(len(flux2)):
             if not (np.isnan(flux2[q])):
             #if not (flux2[q] == flux2[q]):
                 indxs.append(q)

         phase2 = phasetmp[indxs]
         #phase3 = phasetmp3
         flux2 = flux2[indxs]
             
         #JP data
         ax1.plot(phase, flux, 'k-')

         if i == 0:
             pshft = 0.0
             merr = 1.0e6
             for pshift in np.linspace(-0.01, 0.01, 500):
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
         phase2 = phase2 + 0.00448
         
             
         #Savitzky-Golay low-pass filtter
         flux2[-1] = flux2[0]
         #flux3 = savgol_filter(flux2, 15, 5, mode='wrap')
         #flux3 = savgol_filter(flux2, 31, 7, mode='wrap')
         flux3 = savgol_filter(flux2, 41, 11, mode='wrap')
         flux3[0:9] = flux2[0:9]
         flux3[-9:-1] = flux2[-9:-1]
         flux3[-1] = flux2[0]

         flux3 = flux2



         #JN data
         #ax1.plot(phase2, flux2, 'r:')
         ax1.plot(phase2, flux2, 'r--')
         #ax1.plot(phase2, flux2, 'r-', linewidth=0.3)
         

         ax1.set_yticks(ax1.get_yticks()[1:-1])

         #frame for the error panel
         ax2 = subplot(gs[(mfiglim+panelh):(mfiglim+panelh+epanelh), i])
         ax2.minorticks_on()
         ax2.set_xlim(xmin, xmax)
         ax2.set_ylim(eymin, eymax)

         if i == 0:
             ax2.set_ylabel('$\Delta$ %',size=lsize)


         if j == 3:
            ax2.set_xlabel('Phase', size=lsize)
            

         ax2.plot([xmin, xmax], [0.0, 0.0], 'r--', linewidth=0.3)


         #interpolate error
         #fluxi = interp1d(phase, flux, kind='linear')
         #fluxi2 = griddata(phase2, flux2, (phase), method='cubic')
         #fluxi2 = griddata(phase2, flux2, (phase), method='linear')
         
         #fluxi = interp1d(phase, flux, kind='cubic')
         #err = (flux/fluxi2 - 1)*100
                  
         #flux2i = interp1d(phase2, flux2, kind='cubic', fill_value='extrapolate')
         #err = (flux/flux2i(phase) - 1)*100

         #for q in range(len(phase)):
         #    print phase[q], err[q], fluxi2[q], flux[q]
             
         #ax2.plot(phase, err, 'k-', linewidth = 0.4)



         #interpolate error from JP
         #fluxi = griddata(phase, flux, (phase2), method='linear')
         fluxi = griddata(phase, flux, (phase2), method='cubic')
         err = (fluxi/flux3 - 1)*100
         ax2.plot(phase2, err, 'k-', linewidth = 0.4)


         #fluxir = griddata(phase, flux, (phase3), method='cubic')
         #errr = (fluxir/fluxr - 1)*100
         #ax2.plot(phase3, errr, 'g-', linewidth = 0.4)


         #optional errors for range of phase shifts
         for pshift in np.linspace(-0.01, 0.01, 10):
             fluxi2 = griddata(phase2+pshift, flux2, (phase), method='cubic')
             err = (flux/fluxi2 - 1)*100
             #ax2.plot(phase, err, 'b-', linewidth = 0.4)

         
    mfiglim += panelh+epanelh+skiph

    


savefig('fig4.pdf', bbox_inches='tight')

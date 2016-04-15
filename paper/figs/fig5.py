import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

def read_JP2_files(fname):
    da = np.genfromtxt(fname, delimiter="    ", comments='#')
    return da[:,0], da[:,1], da[:,2], da[:,3],da[:,4],da[:,5]

def read_JP_files(fname):
    da = np.genfromtxt(fname, delimiter="   ")
    return da[:,0], da[:,1], da[:,2], da[:,3],da[:,4],da[:,5]

def read_JN_files(fname):
    da = np.genfromtxt(fname, delimiter=",")
    return da[:,0],da[:,1],da[:,2],da[:,3],da[:,4],da[:,5]

def read_Scott_files(fname):

    xx = []
    yy = []
    
    for line in open(fname):
        #print line.rstrip('\n')[5:12], line.rstrip('\n')[26:37]
        #print float(line.rstrip('\n')[5:12]), float(line.rstrip('\n')[26:39])
        xx.append(float(line.rstrip('\n')[5:12]))
        yy.append(float(line.rstrip('\n')[26:39]))
                  
    return xx, yy
                  

    
## Plot
fig = figure(figsize=(9,10), dpi=80)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(400, 3)
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

path3 = "../../out3/waveforms/"

#labels
tsize = 10.0


#nu = '1'
#nu = '400'

fig.text(0.5, 0.92, '$R = 12$ km $\\nu = 400$ Hz $\\rho = 1^{\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.72, '$R = 12$ km $\\nu = 400$ Hz $\\rho = 30^{\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.52, '$R = 15$ km $\\nu = 600$ Hz $\\rho = 1^{\circ}$',  ha='center', va='center', size=tsize)
fig.text(0.5, 0.32, '$R = 15$ km $\\nu = 600$ Hz $\\rho = 30^{\circ}$',  ha='center', va='center', size=tsize)



for j in range(4):

    if j == 0:
        fname = path3 + 'small-1-'
        #fname2 = path3 + 'f400pbbr12m1.6d50i60x1.csv'
        fname2 = path3 + 'f400pbbr12m1.6d50i60x1_+q.csv'
    if j == 1:
        fname = path3 + 'small-30-'
        #fname2 = path3 + 'f400pbbr12m1.6d50i60x30.csv'
        fname2 = path3 + 'f400pbbr12m1.6d50i60x30_+q.csv'
    if j == 2:
        fname = path3 + 'large-1-'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x1.csv'
        fname2 = path3 + 'f600pbbr15m1.6d50i60x1_+q.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x1.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x1_HT+q+w+b.csv'
        #fname2 = path3 + 'r15x1d50i60f600.txt'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x1_Kr.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x1_sHT.csv'
    if j == 3:
        fname = path3 + 'large-30-'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30.csv'
        fname2 = path3 + 'f600pbbr15m1.6d50i60x30_+q.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30_HT.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30_HT-q+w+b.csv' #normal
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30_HT+q+w+b.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30_HT+q+w-b.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30_HT-q-w+b.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30_Kr.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30_exact.csv'
        #fname2 = path3 + 'f600pbbr15m1.6d50i60x30_-Kr.csv'

    #read JP data
    phasea, N2kev = read_Scott_files(fname+'2')
    phaseb, N6kev = read_Scott_files(fname+'6')
    phasec, N12kev = read_Scott_files(fname+'12')
    
    #read JN data
    #phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2 = read_JN_files(fname2) 

    #if j == 1:
    #    phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2 = read_JP2_files(fname2)
    #else:
    phase2, N2kev2, N6kev2, N12kev2, Nbol2, Fbol2 = read_JN_files(fname2)
    phasetmp = phase2

    #bonus file for comparison
    #phase3, N2kev3, N6kev3, N12kev3, Nbol3, Fbol3 = read_JN_files(fname3)
    #phasetmp3 = phase3

    for i in range(3):


         #frame for the main pulse profile fig
         ax1 = subplot(gs[mfiglim:mfiglim+panelh, i])
         ax1.minorticks_on()
         ax1.set_xticklabels([])
         ax1.set_xlim(xmin, xmax)

         formatter = ScalarFormatter(useMathText=True)
         formatter.set_scientific(True)
         formatter.set_powerlimits((0,0))
         ax1.yaxis.set_major_formatter(formatter)
         
         
         if i == 0:
             ax1.set_ylabel('$N$ (2 keV)\n[ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$]',size=lsize)
             phase = phasea
             flux = N2kev
             flux2 = N2kev2
             #flux3 = N2kev3
         elif i == 1:
             ax1.set_ylabel('$N$ (6 keV)',size=lsize)
             phase = phaseb
             flux = N6kev
             flux2 = N6kev2
             #flux3 = N6kev3
         elif i == 2:
             ax1.set_ylabel('$N$ (12 keV)',size=lsize)
             phase = phasec
             flux = N12kev
             flux2 = N12kev2
             #flux3 = N12kev3
         #elif i == 3:
         #    ax1.set_ylabel('Bolometric [ph cm$^{-2}$ s$^{-1}$]',size=lsize)
         #    #flux = Nbol
         #    flux2 = Nbol2

             
         #xxx
         #flux2 = flux
         #phase2 = phase
             
         #Scott data
         ax1.plot(phase, flux, 'k-')

         
         #JN data
         indxs = []
         for q in range(len(flux2)):
             if not (np.isnan(flux2[q])):
                 indxs.append(q)

         phase2 = phasetmp[indxs]
         flux2 = flux2[indxs]
         #phase3 = phasetmp3

         if i == 0:
             pshft = 0.0
             merr = 1.0e6
             for pshift in np.linspace(-0.1, 0.1, 100):
                 fluxi2 = griddata(phase2 + pshift, flux2, (phase), method='cubic', fill_value=0.0)
                 #err = (flux/fluxi2 - 1)*100
                 err = (fluxi2/flux - 1)*100
                 #if err[0] < merr:

                 #print "pshift",pshift," err:"
                 #print err
                 #print

                 serr = 0.0
                 for ijk in range(len(err)):
                     if fluxi2[ijk] != 0:
                         serr += np.abs(err[ijk])

                 
                 #if np.sum(np.abs(err)) < merr: 
                 if serr < merr:
                     #merr = err[0]
                     #merr = np.sum(np.abs(err))
                     merr = serr
                     pshft = pshift
             
             print "min shift:", pshft

         if j == 0:
              phase2 = phase2 
         elif j == 1:
              phase2 = phase2 
         elif j == 2:
              phase2 = phase2 + 0.077 - pshft
         elif j == 3:
             phase2 = phase2 + 0.077 - pshft


         phase2 = phase2 + pshft    
         #phase3 = phase2
         
         ax1.plot(phase2, flux2, 'r--')
         #ax1.plot(phase3, flux3, 'b--')
         
         #frame for the error panel
         ax2 = subplot(gs[(mfiglim+panelh):(mfiglim+panelh+epanelh), i])
         ax2.minorticks_on()
         ax2.set_xlim(xmin, xmax)
         ax2.set_ylim(eymin, eymax)

         if i == 0:
             ax2.set_ylabel('$\Delta$ %',size=lsize)

         #if j != 3:
         #   ax2.set_xticklabels([])

         if j == 3:
            ax2.set_xlabel('Phase', size=lsize)
            
         ax2.plot([xmin, xmax], [0.0, 0.0], 'r--', linewidth=0.3)

         #interpolate error
         fluxi2 = griddata(phase2, flux2, (phase), method='cubic')
         err = (flux/fluxi2 - 1)*100
                  
         ax2.plot(phase, err, 'k-', linewidth = 0.4)


         for pshift in np.linspace(-0.01, 0.01, 10):
             fluxi2 = griddata(phase2+pshift, flux2, (phase), method='cubic')
             err = (fluxi2/flux - 1)*100
             #ax2.plot(phase, err, 'b-', linewidth = 0.4)



         
    mfiglim += panelh+epanelh+skiph

    


savefig('fig5.pdf', bbox_inches='tight')
#savefig('fig2b.pdf', bbox_inches='tight')

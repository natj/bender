import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

def read_csv_files(fname):
    da = np.genfromtxt(fname, delimiter=",")
    return da[:,0], da[:,1]

def read_JN_files(fname):
    da = np.genfromtxt(fname, delimiter=",")
    return da[:,0],da[:,1],da[:,2],da[:,3],da[:,4],da[:,5]


    
## Plot
fig = figure(figsize=(9,4), dpi=80)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(100, 2)
#gs.update(wspace = 0.34)
#gs.update(hspace = 0.4)


lsize = 7.0

xmin = -0.04
xmax = 1.04

eymin = -2.0
eymax = 2.0

panelh = 45
epanelh = 25
skiph = 30


path_Mor = "../../out2/cadeau+morsink/"

#labels
tsize = 10.0

mfiglim = 0

fig.text(0.5, 0.92, '$\\nu = 600$ Hz  blackbody  $\\rho = 1^{\circ}$',  ha='center', va='center', size=tsize)




filesuffix = ['SD', 'SDobl', 'exact']
filelines = ['k-', 'b-', 'r-']


for j in range(1):
    for i in range(2):

        #frame for the main pulse profile fig
         ax1 = subplot(gs[mfiglim:mfiglim+panelh, i])
         ax1.minorticks_on()
         ax1.set_xticklabels([])
         ax1.set_xlim(xmin, xmax)
         ax1.set_ylabel('Flux (arb)',size=lsize)

         #frame for the error panel
         ax2 = subplot(gs[(mfiglim+panelh):(mfiglim+panelh+epanelh), i])
         ax2.minorticks_on()
         ax2.set_xlim(xmin, xmax)
         #ax2.set_ylim(eymin, eymax)
         ax2.set_ylim(-2, 2)
         ax2.set_ylabel('$\Delta$ %',size=lsize)
         ax2.set_xlabel('Phase', size=lsize)
         ax2.plot([xmin, xmax], [0.0, 0.0], 'k', linestyle='dotted', linewidth=0.3)


         for q in range(1):

             q += 2
             #read JN data & Morsink data
             if q == 0:
             
                 if i == 0:
                     fname = path_Mor + 'f600pbbr15m1.4d49i70x1_'+filesuffix[q]+'.csv'
                     fname2 = path_Mor + 'Morsink_f600m14r164i70d49x1_'+filesuffix[q]+'.csv'
                 elif i == 1:
                     fname = path_Mor + 'f600pbbr15m1.4d41i20x1_'+filesuffix[q]+'.csv'
                     fname2 = path_Mor + 'Morsink_f600m14r164i20d41x1_'+filesuffix[q]+'.csv'
             else:
                 if i == 0:
                     fname = path_Mor + 'f600pbbr16m1.4d49i70x1_'+filesuffix[q]+'.csv'
                     fname2 = path_Mor + 'Morsink_f600m14r164i70d49x1_'+filesuffix[q]+'.csv'
                 elif i == 1:
                     fname = path_Mor + 'f600pbbr16m1.4d41i20x1_'+filesuffix[q]+'.csv'
                     fname2 = path_Mor + 'Morsink_f600m14r164i20d41x1_'+filesuffix[q]+'.csv'

                     
             phase, N2kev, N6kev, N12kev, Nbol, Fbol = read_JN_files(fname) 
             flux = Fbol

             

             phase2, flux2 = read_csv_files(fname2)
             fscale = flux.max()

             if q == 0:
                 fscale = flux.max()
             elif q == 1:
                 fscale = 0.00405
             #elif q == 2:
             #    fscale = 0.0041
             #fscale = 0.00498
             flux = flux / fscale
             print 'ratio: ', fscale,' ',filesuffix[q],' ',i


             
             for pshift in np.linspace(-0.01, 0.01, 10):
                 fluxi2 = griddata(phase2+pshift, flux2, (phase), method='cubic')
                 err = (fluxi2/flux - 1)*100
                 ax2.plot(phase, err, 'b-', linewidth = 0.4)


             if q == 0:
                 if i == 0:
                     phase = phase - 0.018
                 elif i == 1:
                     phase = phase - 0.026
             if q == 1:
                 if i == 0:
                     phase = phase + 0.0055
                 elif i == 1:
                     phase = phase + 0.01
             if q == 2:
                 if i == 0:
                     phase = phase + 0.0055
                 elif i == 1:
                     phase = phase + 0.007

                     
             #Morsink data
             if q == 0:
                 ax1.plot(phase2, flux2, 'k-')
             elif q == 1:
                 ax1.plot(phase2, flux2, 'b-')
             elif q == 2:
                 ax1.plot(phase2, flux2, '-', color='darkorange')
                 
             #JN data
             ax1.plot(phase, flux, 'r--')
         
             #print phase2
             
             #interpolate error
             #fluxi2 = interp1d(phase2, flux2, kind='cubic', fill_value=1.0)
             fluxi2 = griddata(phase2, flux2, (phase), method='cubic')
             
             #print 'loop'
             #for k in range(len(phase)):
             #    print phase[k], ' ', fluxi2(phase[k])
             
             #err = (fluxi2(phase)/flux - 1)*100
             err = (fluxi2/flux - 1)*100

             if q == 0:
                 ax2.plot(phase, err, 'k-', linewidth=0.4)
             elif q == 1:
                 ax2.plot(phase, err, 'b-', linewidth=0.4)
             elif q == 2:
                 ax2.plot(phase, err, '-', color='darkorange', linewidth=0.4)
             
             #ax2.plot(phase, err, 'k-', linewidth = 0.4)

             
    mfiglim += panelh+epanelh+skiph
    


savefig('fig3.pdf', bbox_inches='tight')

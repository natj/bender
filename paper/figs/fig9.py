import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

from scipy.signal import savgol_filter



#Read JN files
def read_lineprof(fname):
    da = np.genfromtxt(fname, delimiter=",")

    des = np.diff(da[:,0])[2]
    norm = np.sum(des*da[:,1])
    return da[:,0],da[:,1]/norm


    
## Plot
fig = figure(figsize=(4,3), dpi=80)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 1)
#gs.update(wspace = 0.34)
#gs.update(hspace = 0.4)


lsize = 10.0

xmin = 0.68
xmax = 0.82

#error window limits
eymin = -0.5
eymax = 0.5


#path to files
path_JN = "../../out3/lines/"

#labels size
tsize = 10.0


nu = '700'

#fig.text(0.5, 0.92, '$\\theta_s = 18^{\\circ}$',  ha='center', va='center', size=tsize)
#fig.text(0.5, 0.72, '$\\theta_s = 45^{\\circ}$',  ha='center', va='center', size=tsize)
#fig.text(0.5, 0.52, '$\\theta_s = 90^{\\circ}$',  ha='center', va='center', size=tsize)
#fig.text(0.5, 0.32, 'Hopf $\\theta_s = 45^{\circ}$',  ha='center', va='center', size=tsize)
#fig.text(0.5, 0.12, 'Phase',ha='center', va='center', size=lsize)

ax1 = subplot(gs[0,0])
ax1.minorticks_on()
ax1.set_xlim(xmin, xmax)
#ax1.set_ylim(0.0, 1.05)

ax1.set_ylabel('Normalized flux',size=lsize)
ax1.set_xlabel('Energy $E/E_0$',size=lsize)


#xx1, yy1 = read_lineprof(path_JN+'lineprof_f700pbbr10m1.4i20.csv')
#ax1.plot(xx1, yy1, "k--")

#xx2, yy2 = read_lineprof(path_JN+'lineprof_obl_HTq0_f700pbbr10m1.4i20.csv')
#ax1.plot(xx2, yy2, "k-")

xx3, yy3 = read_lineprof(path_JN+'lineprof_obl_HTq1_f700pbbr10m1.4i20.csv')
ax1.plot(xx3, yy3, "k-", label="$q = -0.268$")

xx4, yy4 = read_lineprof(path_JN+'lineprof_obl_HTq2_f700pbbr10m1.4i20.csv')
ax1.plot(xx4, yy4, "r-", label="$q \\times 2$")

xx5, yy5 = read_lineprof(path_JN+'lineprof_obl_HTq3_f700pbbr10m1.4i20.csv')
ax1.plot(xx5, yy5, "g-", label="$q \\times 3$")

xx6, yy6 = read_lineprof(path_JN+'lineprof_obl_HTq4_f700pbbr10m1.4i20.csv')
ax1.plot(xx6, yy6, "b-", label="$q \\times 4$")

xx7, yy7 = read_lineprof(path_JN+'lineprof_obl_HTq5_f700pbbr10m1.4i20.csv')
ax1.plot(xx7, yy7, "m-", label="$q \\times 5$")

legend = ax1.legend(loc='upper left', shadow=False, labelspacing=0.1)
for label in legend.get_texts():
    label.set_fontsize('x-small')

savefig('fig9.pdf', bbox_inches='tight')

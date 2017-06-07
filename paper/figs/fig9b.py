import numpy as np
import math
from pylab import *

from palettable.wesanderson import Zissou_5 as wsZ
import matplotlib.ticker as mtick

from scipy.interpolate import interp1d
from scipy.interpolate import griddata

from scipy.signal import savgol_filter


from matplotlib import cm
import palettable as pal





#Read JN files
def read_lineprof(fname):
    da = np.genfromtxt(fname, delimiter=",")

    des = np.diff(da[:,0])[2]
    norm = np.sum(des*da[:,1])
    return da[:,0],da[:,1]/norm


#Read JN files
def read_csv(fname):
    da = np.genfromtxt(fname, delimiter=",")

    des = np.diff(da[:,0])[2]
    norm = np.sum(des*da[:,1])
    return da[:,0],da[:,1] #/norm
    
## Plot
fig = figure(figsize=(5,3), dpi=80)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 1)
#gs.update(wspace = 0.34)
#gs.update(hspace = 0.4)


lsize = 10.0

xmin = 0.69
xmax = 0.82


#path to files
path_JN = "../../out/"

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
ax1.set_ylim(0.0, 1.05)

ax1.set_ylabel('Normalized flux',size=lsize)
ax1.set_xlabel('Energy $E/E\'$',size=lsize)


#cmap = pal.cmocean.sequential.Matter_20.mpl_colormap #best so far
cmap = cm.get_cmap('inferno')



files = ['lineprofile_f700_bb_r10_m1.4_i11.csv',
         'lineprofile_f700_bb_r10_m1.4_i10.csv',
         'lineprofile_f700_bb_r10_m1.4_i09.csv',
         'lineprofile_f700_bb_r10_m1.4_i08.csv',
         'lineprofile_f700_bb_r10_m1.4_i07.csv',
         'lineprofile_f700_bb_r10_m1.4_i06.csv',
         'lineprofile_f700_bb_r10_m1.4_i05.csv'
         #'lineprofile_f700_bb_r10_m1.4_i10.csv'
        ]


for i, fname in enumerate(files):
    xx, yy = read_lineprof(path_JN+fname)

    col = cmap( i/float(len(files)) )

    yy /= np.max(yy)
    ax1.plot(xx, yy, "-", linestyle ='solid', color=col)


#cols = ["black",
#        "blue",
#        "red",
#        "magenta"]



path_Bau = "../../out/bau/"

files_Bau = [
        "obl_HT_fit.csv"
        ]


for i, fname in enumerate(files_Bau):
    xx, yy = read_lineprof(path_Bau+fname)

    yy /= np.max(yy)
    ax1.plot(xx, yy, "b-", linestyle ='dashed')




savefig('fig9_b.pdf', bbox_inches='tight')

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

xmin = 0.68
xmax = 0.95

#error window limits
eymin = -0.5
eymax = 0.5


#path to files
path_JN = "../../out3/lines/sweep/"

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
#ax1.set_ylim(0.0, 15)

ax1.set_ylabel('Normalized flux',size=lsize)
ax1.set_xlabel('Energy $E/E\'$',size=lsize)


#xx1, yy1 = read_lineprof(path_JN+'lineprof_f700pbbr10m1.4i20.csv')
#ax1.plot(xx1, yy1, "k--")

#xx2, yy2 = read_lineprof(path_JN+'lineprof_obl_HTq0_f700pbbr10m1.4i20.csv')
#ax1.plot(xx2, yy2, "k-")

#lineprof_obl_HTq3_f700pbbr10m1.4i20.csv
#lineprof_obl_HTq5_f700pbbr10m1.4i20.csv
#lineprof_obl_HTq2_f700pbbr10m1.4i20.csv

files_JN = [
"lineprof_f700pbbr10m1.4i20.csv",
"lineprof_obl_f700pbbr10m1.4i20.csv",
#"lineprof_sph2_HTqfix_f700pbbr10m1.4i20.csv"]
#"lineprof_obl_HTq0_f700pbbr10m1.4i20.csv",
"lineprof_obl_HTq1_f700pbbr10m1.4i20.csv"]
#"lineprof_obl_HTq4_f700pbbr10m1.4i20.csv"]

cols = ["black",
        "blue",
        "red",
        "magenta"]

files_JN = [
'lineprof_obl_HTq1_f700pbbr10m1.4i20.csv',
'lineprof_obl_HTq1_f700pbbr11m1.4i20.csv',
'lineprof_obl_HTq1_f700pbbr12m1.4i20.csv',
'lineprof_obl_HTq1_f700pbbr13m1.4i20.csv',
'lineprof_obl_HTq1_f700pbbr14m1.4i20.csv',
'lineprof_obl_HTq1_f700pbbr15m1.4i20.csv',
'lineprof_obl_HTq1_f700pbbr16m1.4i20.csv']
#'lineprof_obl_HTq1_f700pbbr17m1.4i20.csv']

files_JN2 = [
#'lineprof_obl_HTq1_f700pbbr9m1.1i20.csv ',
'lineprof_obl_HTq1_f700pbbr10m1.1i20.csv',
'lineprof_obl_HTq1_f700pbbr11m1.1i20.csv',
'lineprof_obl_HTq1_f700pbbr12m1.1i20.csv',
'lineprof_obl_HTq1_f700pbbr13m1.1i20.csv',
'lineprof_obl_HTq1_f700pbbr14m1.1i20.csv',
'lineprof_obl_HTq1_f700pbbr15m1.1i20.csv',
'lineprof_obl_HTq1_f700pbbr16m1.1i20.csv']
#'lineprof_obl_HTq1_f700pbbr17m1.1i20.csv']


i = 0
for file_name in files_JN:
    xx, yy = read_lineprof(path_JN+file_name)
    #yy = yy + 3*i
    ax1.plot(xx, yy, color='black', linestyle="solid")
    i += 1



#i = 0
#for file_name in files_JN2:
#    xx, yy = read_lineprof(path_JN+file_name)
#    #yy = yy + 3*i
#    yy = -1.0 * yy
#    ax1.plot(xx, yy, color='red', linestyle="solid")
#    i += 1

ax1.text(0.805, 3, "$R=10$ km", rotation=-85, ha='center', va='center', size=8)
ax1.text(0.832, 3, "$R=11$ km", rotation=-85, ha='center', va='center', size=8)
ax1.text(0.854, 3, "$R=12$ km", rotation=-83, ha='center', va='center', size=8)
ax1.text(0.872, 3, "$R=13$ km", rotation=-82, ha='center', va='center', size=8)
ax1.text(0.887, 3, "$R=14$ km", rotation=-81, ha='center', va='center', size=8)
ax1.text(0.9, 3, "$R=15$ km", rotation=-78, ha='center', va='center', size=8)
ax1.text(0.912, 3, "$R=16$ km", rotation=-77, ha='center', va='center', size=8)


#
#legend = ax1.legend(loc='upper left', shadow=False, labelspacing=0.1)
#for label in legend.get_texts():
#    label.set_fontsize('x-small')

savefig('fig9c.pdf', bbox_inches='tight')

import numpy as np
import matplotlib as mpl
from pylab import *
from matplotlib import cm
from matplotlib.colors import LogNorm







mpl.rcParams['image.cmap'] = 'inferno'
mpl.rc('font', family='serif')
mpl.rc('xtick', labelsize='small')
mpl.rc('ytick', labelsize='small')
gs = GridSpec(1, 3)
gs.update(hspace = 0.3)



#Construct output xy image plane from img object
##################################################
x_span = 11.0
y_span = 11.0

x_bins = 500
y_bins = 500


xs = np.linspace(-x_span, x_span, x_bins)
ys = np.linspace(-y_span, y_span, y_bins)



##################################################
# plot values on image plane
def trans(mat):
    return np.flipud(mat.T)
#return mat
def detrans(mat):
    return np.flipud(mat).T

def clean_image(mat):

    #mask all 0.0 elements and transpose
    mat_masked = np.ma.masked_where(mat == 0, mat) 

    return trans(mat_masked)



#read redshift array
fname = "reds_f600pbbr15m1.4i45.csv"
data = np.genfromtxt(fname, delimiter=',')
redshift = np.reshape(data, (x_bins, y_bins) )
redshift = clean_image(redshift)


##################################################
fname2 = 'reds_f600_bb_r15_m1.4_i45.csv'
data2 = np.genfromtxt(fname2, delimiter=',')
redshift2 = np.reshape(data2, (x_bins, y_bins) )
redshift2 = clean_image(redshift2)






# other settings for imshow
extent=( xs[0], xs[-1], ys[0], xs[-1] )
interpolation = 'nearest'


################################################### 
ax = subplot(gs[0])
ax.minorticks_on()

cax = ax.imshow(redshift, interpolation=interpolation, origin='lower', extent=extent,
        cmap=cm.get_cmap('coolwarm_r'))
ax.contour(redshift, 20, hold='on', colors='w',
        origin='lower', extent=extent)

###################################################
ax = subplot(gs[1])
ax.minorticks_on()

cax = ax.imshow(redshift2, interpolation=interpolation, origin='lower', extent=extent,
        cmap=cm.get_cmap('coolwarm_r'))
ax.contour(redshift2, 20, hold='on', colors='w',
        origin='lower', extent=extent)


ax = subplot(gs[2])
ax.minorticks_on()


###################################################
# relative error
relerr = np.zeros(( x_bins, y_bins))
for i, x in enumerate(xs):
    for j, y in enumerate(ys):

        val1 = redshift[i,j]
        val2 = redshift2[i,j]

        errval = 0.0
        if not(val2 == 0.0):
            errval = np.abs( (val2 - val1)/val2 )
            #errval = np.log10( np.abs((val2 - val1)/val2) )
        relerr[i,j] = errval

relerr = np.ma.masked_where(relerr == 0, relerr) 



#emin = -0.02
#emax =  0.02

print "min :",np.min(relerr)
print "max :",np.max(relerr)

#emin = -3.0
#emax =  1.0


emin = 1.0e-4
emax = 1.0e-1

cax = ax.imshow(relerr,
        interpolation=interpolation, 
        origin='lower', extent=extent,
        cmap=cm.get_cmap('inferno_r'),
        norm=LogNorm(emin, emax)
        #vmin = emin,
        #vmax = emax,
        )



levels = np.linspace(emin, emax, 10)
#levels = np.array( [1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2, 1.0e-1, 5.0e-1, 1.0e0] )
#levels = np.array( [1.0e-2, 2.0e-2 ] )
levels = np.array( [1.0e-3, 5.0e-3] )

ax.contour(relerr, 
        levels, 
        hold='on', 
        linestyle='dashed',
        colors='r',
        origin='lower', 
        extent=extent,
        vmin = emin,
        vmax = emax
        )
colorbar(cax)


show()


savefig('reds.pdf')

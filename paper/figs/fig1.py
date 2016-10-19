import numpy as np
import math
from pylab import *
import matplotlib.patches as patches
from scipy.optimize import minimize_scalar


## Plot
fig = figure(figsize=(6,4), dpi=80)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 1)
#gs.update(wspace = 0.34)
#gs.update(hspace = 0.4)

ax = subplot(gs[0,0])
ax.axis('off')
ax.set_aspect('equal')

#image dimensions
mlim= 1.3
xmin = -1.5
xmax = mlim
ymin = -0.5
ymax = mlim

xlen = xmax-xmin
ylen = ymax-ymin
aspr=ylen/xlen

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)


#variables
#incl = pi/8
incl = pi/2.5
req = 1
u = 0.0
rg = u*req
muc = -rg/(req-rg)
o21 = 0.3

def R(theta):
    return 1.0

def R(theta):
    return 1.0*(1-o21*cos(theta)**2)

def dR(theta):
    return 2.0*o21*sin(theta)*cos(theta)

def fa(theta):
    return dR(theta)/R(theta)

def cosg(theta):
    return (1/sqrt(1-u))/sqrt(1+fa(theta)**2)

def sing(theta):
    return fa(theta)*cosg(theta)

def mu(phi, theta):
    return cos(incl)*cos(theta)+cos(phi)*sin(incl)*sin(theta)

def bimp(phi, theta):
    return sqrt(req*(req+rg+req*mu(phi,theta)-rg*mu(phi,theta))/(1+mu(phi,theta)))

def x(phi, theta):
    return bimp(phi,theta)*R(theta)*(cos(incl)*cos(theta)+cos(phi)*sin(incl)*sin(theta))

def y(phi, theta):
    return bimp(phi,theta)*R(theta)*(sin(theta)*sin(phi))

def z(phi, theta):
    return  bimp(phi,theta)*R(theta)*(cos(theta)*sin(incl)-cos(incl)*cos(phi)*sin(theta))


def xfull(phi,theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (cos(incl)*(a*cos(theta)-c*sin(theta)) + sin(incl)*(c*cos(theta)*cos(phi) + a*cos(phi)*sin(theta) - b*sin(phi)))

def yfull(phi, theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (b*cos(phi) + (c*cos(theta) + a*sin(theta))*sin(phi))

def zfull(phi, theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (sin(incl)*(a*cos(theta)-c*sin(theta)) + cos(incl)*(-cos(phi)*(c*cos(theta) + a*sin(theta)) + b*sin(phi)))



def draw_longitude(ax,
                   phi,
                   start=0,
                   stop=pi/2,
                   rfac=1.0,
                   backside=False,
                   fmt={'color':'k','linestyle':'solid',},
                   ):

    xx = []
    yy = []

    for theta in np.linspace(start, stop, 200):
        if not(backside):
            if mu(phi,theta) >= muc:
                xx.append(y(phi,theta)*rfac)
                yy.append(z(phi,theta)*rfac)
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []
        else:
            if mu(phi,theta) <= muc:
                xx.append(y(phi,theta)*rfac)
                yy.append(z(phi,theta)*rfac)
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []

    lines, = ax.plot(xx,yy)
    lines.set(**fmt)

    return ax


def draw_latitude(ax,
                   theta,
                   start=-pi,
                   stop=pi,
                   rfac=1.0,
                   backside=False,
                   fmt={'color':'k','linestyle':'solid',},
                   ):

    xx = []
    yy = []

    for phi in np.linspace(start, stop, 200):
        if not(backside):
            if mu(phi,theta) >= muc:
                xx.append(y(phi,theta)*rfac)
                yy.append(z(phi,theta)*rfac)
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []
        else:
            if mu(phi,theta) <= muc:
                xx.append(y(phi,theta)*rfac)
                yy.append(z(phi,theta)*rfac)
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []

    lines, = ax.plot(xx, yy)
    lines.set(**fmt)

    return ax

def draw_outline(ax,
                 fmt={'color':'k','linestyle':'solid',},
                 ):

    def rootf(theta1, phi1):
        return abs(mu(phi1, theta1))-muc

    xx = []
    yy = []

    #print "outline:"
    for phi in np.linspace(pi/2, 3*pi/2, 100):
        #print "phi:", phi
        thetamin = pi/2
        thetamax = 0.0

        res = minimize_scalar(rootf, bounds=(thetamin, thetamax), args=(phi,))
        #print "minim:", res.x, " ", res.fun

        theta2 = res.x

        #ax.plot([y(phi, theta2)], [z(phi, theta2)], "b.")
        xx.append( y(phi, theta2))
        yy.append( z(phi, theta2))

    ax.plot(xx, yy, "r-")

    return ax


def draw_radial(ax, phi, theta,
                fmt={'color':'k','linestyle':'solid',},
                rfac=1.0):
    xx = [0.0]
    yy = [0.0]

    xx.append(y(phi, theta)*rfac)
    yy.append(z(phi, theta)*rfac)

    lines, = ax.plot(xx, yy)
    lines.set(**fmt)

    return ax


def draw_axis(ax, phi, theta,
              startp=(0.0, 0.0),
              fmt= {'color':'k', 
              'linestyle':'solid', 
              'lw':1.5, 
              'head_width': 0.04, 
              'head_length': 0.08,},
              rfac=1.0):

    ax.add_patch(patches.FancyArrow(
              startp[0],startp[1],
              y(phi, theta)*rfac, z(phi, theta)*rfac,
              **fmt
              ))

    return ax


#normal vector
def draw_normal(ax, phi, theta,
                fmt={'color':'black', 'linestyle':'solid', 'lw':1.8, 'head_width': 0.02, 'head_length': 0.05,},
                alen=0.3
                ):

    spot_cosg = cosg(theta)
    spot_sing = sing(theta)

    normal_y1 = y(phi, theta)*1.0 
    normal_z1 = z(phi, theta)*1.0,
    normal_y2 = yfull(phi, theta, spot_cosg, 0.0, -spot_sing)*alen
    normal_z2 = zfull(phi, theta, spot_cosg, 0.0, -spot_sing)*alen

    ax.add_patch(patches.FancyArrow(
                normal_y1,      normal_z1,
                normal_y2,      normal_z2,
                **fmt
                ))
    
    return ax


def draw_observer(ax, xx, yy, angle, sangle=0.3, alen=0.2, 
                  fmt={'color':'k','linestyle':'solid','lw':1.0}):

    x1 = xx + alen*cos(angle-sangle/2)
    y1 = yy + alen*sin(angle-sangle/2)

    x2 = xx + alen*cos(angle+sangle/2)
    y2 = yy + alen*sin(angle+sangle/2)

    ax.plot([xx,x1],[yy,y1],**fmt)
    ax.plot([xx,x2],[yy,y2],**fmt)

    xsegm = []
    ysegm = []
    for theta in np.linspace(angle-sangle/2-0.1, angle+sangle/2+0.1, 10):
        xsegm.append(xx + alen*0.7*cos(theta))
        ysegm.append(yy + alen*0.7*sin(theta))
    ax.plot(xsegm, ysegm, **fmt)

    return ax

def draw_spot(ax, sphi, stheta, rho,
              fmt={'facecolor':'k', 'alpha':0.3,}):

    #TODO check visibility using mu

    xx = []
    yy = []
    for chi in np.linspace(0, 2*pi, 30):

        #transform from spot centric-coordinate system into the standard spherical coord
        xs = cos(rho)*cos(sphi)*sin(stheta)+sin(rho)*(cos(stheta)*cos(sphi)*cos(chi)-sin(sphi)*sin(chi))
        ys = cos(rho)*sin(stheta)*sin(sphi)+sin(rho)*(cos(stheta)*cos(chi)*sin(sphi)+cos(sphi)*sin(chi))
        zs = cos(stheta)*cos(rho)-cos(chi)*sin(stheta)*sin(rho)

        #transform to phi, theta
        phi = arctan2(ys,xs)
        theta = arccos(zs/sqrt(xs**2 + ys**2 + zs**2))

        #ax.plot([y(phi,theta)], [z(phi,theta)], "b.")

        xx.append(y(phi, theta))
        yy.append(z(phi, theta))

    #ax.plot(xx, yy, "b-")
    ax.add_patch(Polygon(zip(xx,yy),
                      closed=True, 
                      #facecolor='black',alpha=0.5))
                      **fmt))

    return ax


#--------------------------------------------------------------------------------



#coordinates of spot and observer
spot_theta=pi/5
spot_phi=0.67

#obs_theta=incl
obs_theta=pi/2.7
obs_phi=-pi/1.5

#default front and back line styles
fmt={'color':'k','linestyle':'solid',}
fmtb = {'color':'k','linestyle':'dotted',}



#-----------------------------------------------------
#draw xyz axis
ax = draw_axis(ax, obs_phi,      pi/2, rfac=1.4)
ax = draw_axis(ax, obs_phi+pi/2, pi/2, rfac=1.4)
ax = draw_axis(ax, 0.0,          0.0 , rfac=1.8)

#rotation direction, i.e. curly arrow
phistart=-pi/2-0.8
phistop=pi/2+0.5
wsize=0.15
wheight=1.2
ax=draw_latitude(ax, wsize, fmt=fmt, start=phistart, stop=phistop, rfac=wheight)
yy1=y(phistop,wsize)*wheight
zz1=z(phistop,wsize)*wheight

yy2 = -0.01
zz2 = 0.002

fmta= {'color':'black', 'linestyle':'solid', 'lw':1.0, 'head_width': 0.02, 'head_length': 0.04,}
ax.add_patch(patches.FancyArrow(
            yy1, zz1,
            yy2, zz2,
            **fmta
            ))



#-----------------------------------------------------
#borders
fmtO={'color':'red','linestyle':'solid',}
#ax = draw_outline(ax, fmtO)

ax=draw_longitude(ax, pi/2, fmt=fmt)
ax=draw_longitude(ax, -pi/2, fmt=fmt)

#ax=draw_longitude(ax, -0.1, fmt=fmt)
#ax=draw_longitude(ax, -pi-0.1, fmt=fmt)
#ax=draw_longitude(ax, -pi-0.1, fmt=fmtb, backside=True)

ax=draw_latitude(ax, pi/2, fmt=fmt)
ax=draw_latitude(ax, pi/2, backside=True, fmt=fmtb)


#fmt2 = {'color':'k','linestyle':'solid','lw':0.8,}
#fmt2b = {'color':'k','linestyle':'dotted','lw':0.8,}
#for thetaa in np.linspace(0.0, pi/2, 10):
#    ax = draw_latitude(ax, thetaa, fmt=fmt2)
#    ax = draw_latitude(ax, thetaa, fmt=fmt2b, backside=True)
#                       
#for phii in np.linspace(0.0, 2*pi, 30):
#    ax = draw_longitude(ax, phii, fmt=fmt2,
#                        start=0.0,
#                        stop=pi/2
#                        )


#-----------------------------------------------------
#spot position
ax=draw_longitude(ax, spot_phi, fmt=fmt)
ax=draw_radial(ax, spot_phi, spot_theta, fmt)
#ax = draw_axis(ax, spot_phi, spot_theta, rfac=1.2)

fmta= {'color':'black', 'linestyle':'solid', 'lw':1.4, 'head_width': 0.02, 'head_length': 0.05,}
alen = 0.6
ax.add_patch(patches.FancyArrow(
            y(spot_phi, spot_theta)*1.0, z(spot_phi, spot_theta)*1.0,
            y(spot_phi, spot_theta)*alen, z(spot_phi, spot_theta)*alen,
            **fmta
            ))


fmta= {'color':'black', 'linestyle':'solid', 'lw':1.4, 'head_width': 0.02, 'head_length': 0.05,}
ax = draw_normal(ax, spot_phi, spot_theta, fmt=fmta, alen=0.4)

#for spot_theta in np.linspace(0.1, pi/2, 10):
#    ax = draw_normal(ax, spot_phi, spot_theta, alen=0.1)


ax=draw_longitude(ax, spot_phi, fmt=fmt, start=0, stop=spot_theta, rfac=0.17)
ax=draw_longitude(ax, spot_phi, fmt=fmt, start=0, stop=spot_theta, rfac=0.17, backside=True) #theta angle
#add(p, PlotLabel(.53, .58, "<i>\\theta</i>"))

#help axis
ax = draw_radial(ax, spot_phi, pi/2, fmt)


#circular spot
rho = 0.1
ax = draw_spot(ax, spot_phi, spot_theta, rho)

#ax=draw_longitude(ax, spot_phi-rho, fmt=fmt)
#ax=draw_longitude(ax, spot_phi+rho, fmt=fmt)
#ax=draw_radial(ax, spot_phi, spot_theta+rho, fmt)
#ax=draw_radial(ax, spot_phi, spot_theta-rho, fmt)


#-----------------------------------------------------
#observer
fmti = {'color':'k','linestyle':'dashed',}

ax=draw_radial(ax, obs_phi, obs_theta, fmt=fmti, rfac=1.65)
ax=draw_radial(ax, obs_phi, pi/2, fmt=fmti)

ax=draw_longitude(ax, obs_phi, fmt=fmt, start=0, stop=obs_theta, rfac=0.1) #incl angle
ax=draw_longitude(ax, obs_phi, fmt=fmt, start=0, stop=obs_theta, rfac=0.1, backside=True) #incl angle


ax=draw_latitude(ax, pi/2, fmt=fmt, start=spot_phi, stop=obs_phi, rfac=0.15) #phi angle
ax=draw_latitude(ax, pi/2, fmt=fmt, start=spot_phi, stop=obs_phi, rfac=0.15, backside=True) #phi angle

obspx = y(obs_phi, obs_theta)*1.8
obspy = z(obs_phi, obs_theta)*1.8
#ax = draw_observer(ax, -1.5, 0.8, angle=-0.45, sangle=0.4, alen=0.15)
ax = draw_observer(ax, obspx, obspy, angle=obs_theta-pi/2-0.16, sangle=0.4, alen=0.15)


#observer axis
fmt_obs= {'color':'k', 
        'linestyle':'solid', 
        'lw':1.5, 
        'head_width': 0.02, 
        'head_length': 0.04,}

obspx = y(obs_phi, obs_theta)*1.25
obspy = z(obs_phi, obs_theta)*1.25
ax = draw_axis(ax, obs_phi,      obs_theta,      rfac=0.2, startp=(obspx, obspy),fmt=fmt_obs)
ax = draw_axis(ax, obs_phi+pi/2, pi/2,           rfac=0.2, startp=(obspx, obspy),fmt=fmt_obs)
ax = draw_axis(ax, obs_phi,      obs_theta-pi/2, rfac=0.2, startp=(obspx, obspy),fmt=fmt_obs)


#texts
#-----------------------------------------------------
lsize = 15.0
ax.text(0.0, -0.15, "$\\phi$", va='center', ha='center', size=lsize)
ax.text(0.06, 0.2, "$\\theta$", va='center', ha='center', size=lsize)
ax.text(-0.05, 0.15, "$i$", va='center', ha='center', size=lsize)

ax.text(-0.6, -0.4, "$x$", va='center', ha='center', size=lsize)
ax.text(0.06, 1.1, "$y$", va='center', ha='center', size=lsize)
ax.text(-1.3, 0.1, "$z$", va='center', ha='center', size=lsize)
#ax.text(-0.1, 0.9, "$\\hat{\\Omega}$", va='center', ha='center', size=lsize)

ax.text(0.55, 0.96, "$\\hat{r}$", va='center', ha='center', size=lsize)
ax.text(0.38, 0.92, "$\\hat{n}$", va='center', ha='center', size=lsize)

ax.text(-1.15, 0.5, "$\\hat{x}$", va='center', ha='center', size=lsize)
ax.text(-0.83, 0.88,  "$\\hat{y}$", va='center', ha='center', size=lsize)
#ax.text(-1.2, 0.8,  "$\\hat{z}$", va='center', ha='center', size=lsize)

savefig('fig1.png', bbox_inches='tight')

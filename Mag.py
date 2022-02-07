import simpegPF as MAG
import simpegCoordUtils as Utils

from scipy.constants import mu_0
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from IPython.html.widgets import *
# import ipywidgets as widgets
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# monFile = "data2015/StudentData2015_Monday.csv"
# monData = pd.DataFrame(pd.read_csv(filename, header = 0))

# filename = "data2014/HZrebarProfile.csv"
# data = pd.DataFrame(pd.read_csv(filename, header = 0))
# loc = data["Distance"].values

diameter = 1.4e-2
length = 3.
xlim = np.r_[5., 25.]
rx_h = 1.9

Bincd = 70.205
Bdecd = 16.63
Bigrfd = 54450

# Depth of burial: Monday was 35cm. I believe that Wednesday was ~45cm


class definePrism(object):
    """
        Define a prism and its attributes

        Prism geometry:
            - dx, dy, dz: width, length and height of prism
            - depth : depth to top of prism
            - susc : susceptibility of prism
            - x0, y0 : center of prism in horizontal plane
            - pinc, pdec : inclination and declination of prism
    """

    x0, y0, z0, dx, dy, dz = 0., 0., 0., 1., 1., 1.
    pinc, pdec = 0., 0.


    # Define the nodes of the prism
    @property
    def xn(self):
        xn = np.asarray([-self.dx/2. + self.x0, self.dx/2. + self.x0])

        return xn

    @property
    def yn(self):
        yn = np.asarray([-self.dy/2. + self.y0, self.dy/2. + self.y0])

        return yn

    @property
    def zn(self):
        zn = np.asarray([-self.dz + self.z0, self.z0])

        return zn

    @property
    def xc(self):
        xc = (self.xn[0] + self.xn[1]) / 2.

        return xc

    @property
    def yc(self):
        yc = (self.yn[0] + self.yn[1]) / 2.

        return yc

    @property
    def zc(self):
        zc = (self.zn[0] + self.zn[1]) / 2.

        return zc


def plotProfile(prob2D, x0, data, Binc, Bdec, Bigrf, susc, Q, rinc, rdec):
    if data is 'MonSt':
        filename = "data/Lab1_monday_TA.csv"
    elif data is 'WedSt':
        filename = "data/Lab1_Wednesday_student.csv"
    elif data is 'WedTA':
        filename = "data/Lab1_Wednesday_TA.csv"

    dat = pd.DataFrame(pd.read_csv(filename, header = 0))
    tf  = dat["MAG_MEAN"].values
    std = dat["STDEV"].values
    loc = dat["station"].values
    #teams = dat["Team"].values

    tfa = tf - Bigrf

    p = prob2D.prism

    # nx, ny = 100, 1
    # shape = (nx, ny)

    dx = x0 - prob2D.survey.xylim

    prob2D.survey.profile
    if prob2D.survey.profile == "EW":
        x1, x2, y1, y2 = prob2D.survey.xr[0]-dx, prob2D.survey.xr[-1]-dx, 0., 0.
    elif prob2D.survey.profile == "NS":
        x1, x2, y1, y2 = 0., 0., prob2D.survey.yr[0]-dx, prob2D.survey.yr[-1]-dx
    elif prob2D.survey.profile == "45N":
        x1, x2, y1, y2 = prob2D.survey.xr[0], prob2D.survey.xr[-1], prob2D.survey.yr[0], prob2D.survey.yr[-1]


    x, y = linefun(x1, x2, y1, y2, 100)
    xyz_line = np.c_[x, y, np.ones_like(x)*prob2D.survey.rx_h]

    distance = np.sqrt((x-x1)**2.+(y-y1)**2.)

    xlim = [0,distance[-1]]
    prob1D = MAG.problem()
    srvy1D = MAG.survey()
    srvy1D._rxLoc = xyz_line

    prob1D.prism = p
    prob1D.survey = srvy1D

    prob1D.Bdec, prob1D.Binc, prob1D.Bigrf = Bdec, Binc, Bigrf
    prob1D.Q, prob1D.rinc, prob1D.rdec = Q, rinc, rdec
    prob1D.uType, prob1D.mType = 'tf', 'total'
    prob1D.susc = susc

    # Compute fields from prism
    magi, magr = prob1D.fields()

    #out_linei, out_liner = getField(p, xyz_line, comp, 'total')
    #out_linei = getField(p, xyz_line, comp,'induced')
    #out_liner = getField(p, xyz_line, comp,'remanent')

    # distance = np.sqrt((x-x1)**2.+(y-y1)**2.)

    f = plt.figure(figsize = (10, 5))
    gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    ax1.plot(x0, p.z0, 'ko')
    ax1.text(x0+0.5, p.z0, 'Rebar', color='k')
    ax1.text(xlim[0]+1.,-1.2, 'Magnetometer height (1.9 m)', color='b')
    ax1.plot(xlim, np.r_[rx_h, rx_h], 'b--')

    # magi,magr = getField(p, rxLoc, 'bz', 'total')

    ax1.plot(xlim, np.r_[0., 0.], 'k--')
    ax1.set_xlim(xlim)
    ax1.set_ylim(-2.5, 2.5)

    ax0.scatter(loc,tfa)
    ax0.errorbar(loc,tfa,yerr=std,linestyle = "None",color="k")
    ax0.set_xlim(xlim)
    ax0.grid(which="both")

    ax0.plot(distance, magi, 'b', label='induced')
    ax0.plot(distance, magr, 'r', label='remnant')
    ax0.plot(distance, magi+magr, 'k', label='total')
    ax0.legend(loc=2,fontsize=10)
    # ax[1].plot(loc-8, magnT[::-1], )

    ax1.set_xlabel("Northing (m)")
    ax1.set_ylabel("Depth (m)")

    ax0.set_ylabel("Total field anomaly (nT)")

    ax0.grid(True)
    ax1.grid(True)

    if prob2D.survey.profile == "EW":
        ax1.set_xlabel("Easting (m)")
        ax0.set_xlabel("Easting (m)")

    elif prob2D.survey.profile == "NS":
        ax1.set_xlabel("Northing (m)")
        ax0.set_xlabel("Northing (m)")

    elif prob2D.survey.profile == "45N":
        ax1.set_xlabel("Distance SW-NE (m)")
        ax0.set_xlabel("Distance SW-NE (m)")

            # ax1.invert_yaxis()

    plt.tight_layout()
    plt.show()

    return True


def plotObj3D(p, rx_h, elev, azim, npts2D, xylim,
              profile=None, x0=15., y0=0., fig=None, axs=None, plotSurvey=True,cstring='k'):

    # define the survey area
    surveyArea = (-xylim, xylim, -xylim, xylim)
    shape = (npts2D, npts2D)

    xr = np.linspace(-xylim, xylim, shape[0])
    yr = np.linspace(-xylim, xylim, shape[1])
    X, Y = np.meshgrid(xr, yr)
    Z = np.ones(np.shape(X))*rx_h

    rxLoc = np.c_[Utils.mkvc(X), Utils.mkvc(Y), Utils.mkvc(Z)]

    depth = p.z0
    x1, x2 = p.xn[0]-p.xc, p.xn[1]-p.xc
    y1, y2 = p.yn[0]-p.yc, p.yn[1]-p.yc
    z1, z2 = p.zn[0]-p.zc, p.zn[1]-p.zc

    pinc, pdec = p.pinc, p.pdec

    if fig is None:
        fig = plt.figure(figsize=(7, 7))

    if axs is None:
        axs = fig.add_subplot(111, projection='3d')
    plt.rcParams.update({'font.size': 13})

    axs.set_xlim3d(surveyArea[:2])
    axs.set_ylim3d(surveyArea[2:])
#     axs.set_zlim3d(depth+np.array(surveyArea[:2]))
    axs.set_zlim3d(-surveyArea[-1]*rx_h, 3)

    # Create a rectangular prism, rotate and plot
    block_xyz = np.asarray([[x1, x1, x2, x2, x1, x1, x2, x2],
                           [y1, y2, y2, y1, y1, y2, y2, y1],
                           [z1, z1, z1, z1, z2, z2, z2, z2]])

    # rot = Utils.mkvc(Utils.dipazm_2_xyz(pinc, pdec))

    # xyz = Utils.rotatePointsFromNormals(block_xyz.T, np.r_[0., 1., 0.], rot,
    #                                     np.r_[p.xc, p.yc, p.zc])

    R = Utils.rotationMatrix(pinc, pdec)

    xyz = R.dot(block_xyz).T

#    print xyz
    # Face 1
    coordsforthis=[]

   # pltvec=[0,1,2,3,4]

   
   
   # for ii in pltvec:
    #    coordsforthis.append([xyz[ii, 0]+p.xc,xyz[ii, 1]+p.yc,xyz[ii, 2]+p.zc])

 #   coordsforthis.append([(xyz[:4,0],xyz[:4,1],xyz[:4,2])])
    xyz2=np.copy(xyz)
    xyz2[:,0]=xyz2[:,0]+p.xc
    xyz2[:,1]=xyz2[:,1]+p.yc
    xyz2[:,2]=xyz2[:,2]+p.zc
    coordsforthis.append([(xyz2[0]),(xyz2[1]),(xyz2[2]),(xyz2[3])])
    axs.add_collection3d(Poly3DCollection(coordsforthis, facecolors=cstring))

    # Face 2
    axs.add_collection3d(Poly3DCollection(coordsforthis, facecolors=cstring))

    coordsforthis2=[]
    coordsforthis2.append([(xyz2[0]),(xyz2[1]),(xyz2[5]),(xyz2[4])])
    # Face 3
    axs.add_collection3d(Poly3DCollection(coordsforthis2, facecolors=cstring))

    coordsforthis3=[]
    coordsforthis3.append([(xyz2[3]),(xyz2[2]),(xyz2[6]),(xyz2[7])])
    
    # Face 4
    axs.add_collection3d(Poly3DCollection(coordsforthis3, facecolors=cstring))


    coordsforthis4=[]
    coordsforthis4.append([(xyz2[0]),(xyz2[4]),(xyz2[7]),(xyz2[3])])
    # Face 5
    axs.add_collection3d(Poly3DCollection(coordsforthis4, facecolors=cstring))

    
    coordsforthis5=[]
    coordsforthis5.append([(xyz2[1]),(xyz2[5]),(xyz2[6]),(xyz2[2])])
    # Face 6
    axs.add_collection3d(Poly3DCollection(coordsforthis5, facecolors=cstring))

    axs.set_xlabel('Easting (X; m)')
    axs.set_ylabel('Northing (Y; m)')
    axs.set_zlabel('Depth (Z; m)')
#     axs.invert_zaxis()
#     axs.invert_yaxis()

    if plotSurvey:
        axs.plot(rxLoc[:, 0], rxLoc[:, 1], rxLoc[:, 2], '.g', alpha=0.5)

    if profile == "EW":
        axs.plot(np.r_[surveyArea[:2]], np.r_[0., 0.], np.r_[rx_h, rx_h], 'r-')
    elif profile == "NS":
        axs.plot(np.r_[0., 0.], np.r_[surveyArea[2:]], np.r_[rx_h, rx_h], 'r-')
    elif profile == "45N":
        axs.plot(np.r_[surveyArea[:2]], np.r_[surveyArea[2:]], np.r_[rx_h, rx_h], 'r-')
        # axs.plot(np.r_[surveyArea[2:]], np.r_[0., 0.], np.r_[rx_h, rx_h], 'r-')

    axs.view_init(elev, azim)
 
    if fig is None:
        plt.show()

    return True


def linefun(x1, x2, y1, y2, nx, tol=1e-3):
    dx = x2-x1
    dy = y2-y1

    if np.abs(dx) < tol:
        y = np.linspace(y1, y2, nx)
        x = np.ones_like(y)*x1
    elif np.abs(dy) < tol:
        x = np.linspace(x1, x2, nx)
        y = np.ones_like(x)*y1
    else:
        x = np.linspace(x1, x2, nx)
        slope = (y2-y1)/(x2-x1)
        y = slope*(x-x1)+y1
    return x, y


def plogMagSurvey2D(prob2D, susc, Einc, Edec, Bigrf, comp, irt,  Q, rinc, rdec, fig=None, axs1=None, axs2=None):

    import matplotlib.gridspec as gridspec

    # The MAG problem created is stored in result[1]
    # prob2D = Box.result[1]

    if fig is None:
        fig = plt.figure(figsize=(18*1.5,3.4*1.5))

        plt.rcParams.update({'font.size': 14})
        gs1 = gridspec.GridSpec(2, 7)
        gs1.update(left=0.05, right=0.48, wspace=0.05)

    if axs1 is None:
        axs1 = plt.subplot(gs1[:2, :3])

    if axs2 is None:
        axs2 = plt.subplot(gs1[0, 4:])

    axs1.axis("equal")

    prob2D.Bdec, prob2D.Binc, prob2D.Bigrf = Edec, Einc, Bigrf
    prob2D.Q, prob2D.rinc, prob2D.rdec = Q, rinc, rdec
    prob2D.uType, prob2D.mType = comp, 'total'
    prob2D.susc = susc

    # Compute fields from prism
    b_ind, b_rem = prob2D.fields()

    if irt == 'total':
        out = b_ind + b_rem

    elif irt == 'induced':
        out = b_ind

    else:
        out = b_rem

    X, Y = np.meshgrid(prob2D.survey.xr, prob2D.survey.yr)

    dat = axs1.contourf(X,Y, np.reshape(out, (X.shape)).T, 25)
    cb = plt.colorbar(dat, ax=axs1, ticks=np.linspace(out.min(), out.max(), 5))
    cb.set_label("nT")

    axs1.plot(X, Y, '.k')

    # Compute fields on the line by creating a similar mag problem
    prob2D.survey.profile
    if prob2D.survey.profile == "EW":
        x1, x2, y1, y2 = prob2D.survey.xr[0], prob2D.survey.xr[-1], 0., 0.
    elif prob2D.survey.profile == "NS":
        x1, x2, y1, y2 = 0., 0., prob2D.survey.yr[0], prob2D.survey.yr[-1]
    elif prob2D.survey.profile == "45N":
        x1, x2, y1, y2 = prob2D.survey.xr[0], prob2D.survey.xr[-1], prob2D.survey.yr[0], prob2D.survey.yr[-1]

    x, y = linefun(x1, x2, y1, y2, prob2D.survey.npts2D)
    xyz_line = np.c_[x, y, np.ones_like(x)*prob2D.survey.rx_h]
    # Create problem
    prob1D = MAG.problem()
    srvy1D = MAG.survey()
    srvy1D._rxLoc = xyz_line

    prob1D.prism = prob2D.prism
    prob1D.survey = srvy1D

    prob1D.Bdec, prob1D.Binc, prob1D.Bigrf = Edec, Einc, Bigrf
    prob1D.Q, prob1D.rinc, prob1D.rdec = Q, rinc, rdec
    prob1D.uType, prob1D.mType = comp, 'total'
    prob1D.susc = susc

    # Compute fields from prism
    out_linei, out_liner = prob1D.fields()

    #out_linei, out_liner = getField(p, xyz_line, comp, 'total')
    #out_linei = getField(p, xyz_line, comp,'induced')
    #out_liner = getField(p, xyz_line, comp,'remanent')

    out_linet = out_linei+out_liner

    distance = np.sqrt((x-x1)**2.+(y-y1)**2.)


    axs1.plot(x,y, 'w.', ms=3)

    axs1.text(x[0], y[0], 'A', fontsize=16, color='w')
    axs1.text(x[-1], y[-1], 'B', fontsize=16,
             color='w', horizontalalignment='right')

    axs1.set_xlabel('Easting (X; m)')
    axs1.set_ylabel('Northing (Y; m)')
    axs1.set_xlim(X.min(), X.max())
    axs1.set_ylim(Y.min(), Y.max())
    axs1.set_title(irt+' '+comp)

    axs2.plot(distance, out_linei, 'b.-')
    axs2.plot(distance, out_liner, 'r.-')
    axs2.plot(distance, out_linet, 'k.-')
    axs2.set_xlim(distance.min(), distance.max())

    axs2.set_xlabel("Distance (m)")
    axs2.set_ylabel("Magnetic field (nT)")

    axs2.text(distance.min(), out_linei.max()*0.8, 'A', fontsize = 16)
    axs2.text(distance.max()*0.97, out_linei.max()*0.8, 'B', fontsize = 16)
    axs2.legend(("induced", "remanent", "total"), bbox_to_anchor=(0.5, -0.3))
    axs2.grid(True)
    plt.show()

    return True


def fitline(Box):

    def profiledata(data, Binc, Bdec, Bigrf, x0, depth, susc, Q, rinc, rdec, update):

        prob = Box.result[1]
        prob.prism.z0 = -depth

        return plotProfile(prob, x0, data, Binc, Bdec, Bigrf, susc, Q, rinc, rdec)

    Q = widgets.interactive(profiledata, data=widgets.ToggleButtons(options=['MonSt','WedTA','WedSt']),\
             Binc=widgets.FloatSlider(min=-90.,max=90,step=5,value=90,continuous_update=False),\
             Bdec=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False),\
             Bigrf=widgets.FloatSlider(min=54000.,max=55000,step=10,value=54500,continuous_update=False),\
             x0=widgets.FloatSlider(min=5., max=25., step=0.1, value=15.), \
             depth=widgets.FloatSlider(min=0.,max=2.,step=0.05,value=0.5), \
             susc=widgets.FloatSlider(min=0., max=800.,step=5., value=1.),\
             Q=widgets.FloatSlider(min=0., max=10.,step=0.1, value=0.),\
             rinc=widgets.FloatSlider(min=-180., max=180.,step=1., value=0.),\
             rdec=widgets.FloatSlider(min=-180., max=180.,step=1., value=0.),
             update=widgets.ToggleButton(description='Refresh', value=False) \
             )
    return Q


def ViewMagSurvey2DInd(Box):


    def MagSurvey2DInd(susc, Einc, Edec, Bigrf, comp, irt, Q, rinc, rdec, update):

        # Get the line extent from the 2D survey for now
        prob = Box.result[1]

        return plogMagSurvey2D(prob, susc, Einc, Edec, Bigrf, comp, irt, Q, rinc, rdec)

    out = widgets.interactive (MagSurvey2DInd
                    ,susc=widgets.FloatSlider(min=0,max=200,step=0.1,value=0.1,continuous_update=False) \
                    ,Einc=widgets.FloatSlider(min=-90.,max=90,step=5,value=90,continuous_update=False) \
                    ,Edec=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False) \
                    ,Bigrf=widgets.FloatSlider(min=53000.,max=55000,step=10,value=54500,continuous_update=False) \
                    ,comp=widgets.ToggleButtons(options=['tf','bx','by','bz'])
                    ,irt=widgets.ToggleButtons(options=['induced','remanent', 'total'])
                    ,Q=widgets.FloatSlider(min=0.,max=10,step=1,value=0,continuous_update=False) \
                    ,rinc=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False) \
                    ,rdec=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False) \
                    ,update=widgets.ToggleButton(description='Refresh', value=False) \
                    )
    return out


def Prism(dx, dy, dz, depth, pinc, pdec, npts2D, xylim, rx_h, profile, View_elev, View_azim):
    #p = definePrism(dx, dy, dz, depth,pinc=pinc, pdec=pdec, susc = 1., Einc=90., Edec=0., Bigrf=1e-6)
    p = definePrism()
    p.dx, p.dy, p.dz, p.z0 = dx, dy, dz, -depth
    p.pinc, p.pdec = pinc, pdec

    srvy = MAG.survey()
    srvy.rx_h, srvy.npts2D, srvy.xylim, srvy.profile = rx_h, npts2D, xylim, profile

    # Create problem
    prob = MAG.problem()
    prob.prism = p
    prob.survey = srvy

    return plotObj3D(p, rx_h, View_elev, View_azim, npts2D, xylim, profile=profile), prob

def PrismAndData(dx, dy, dz, depth, pinc, pdec, npts2D, xylim, rx_h, profile, View_elev, View_azim,susc,Binc,Bdec,Bigrf):
    #p = definePrism(dx, dy, dz, depth,pinc=pinc, pdec=pdec, susc = 1., Einc=90., Edec=0., Bigrf=1e-6)
    p = definePrism()
    p.dx, p.dy, p.dz, p.z0 = dx, dy, dz, -depth
    p.pinc, p.pdec = pinc, pdec

    srvy = MAG.survey()
    srvy.rx_h, srvy.npts2D, srvy.xylim, srvy.profile = rx_h, npts2D, xylim, profile

    # Create problem
    prob = MAG.problem()
    prob.prism = p
    prob.survey = srvy

    fig=plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1]) 
    ax1=fig.add_subplot(gs[0],projection='3d')
    ax2=fig.add_subplot(gs[1])
    ax3=0

    prismplt=plotObj3D(p, rx_h, View_elev, View_azim, npts2D, xylim, profile=profile,fig=fig,axs=ax1)

    PlotPipeField(prob,Binc,Bdec,Bigrf,susc,fig=fig,axs=[ax2,ax3],dolegend=False)


    
def linefun(x1, x2, y1, y2, nx, tol=1e-3):
    dx = x2-x1
    dy = y2-y1

    if np.abs(dx) < tol:
        y = np.linspace(y1, y2, nx)
        x = np.ones_like(y)*x1
    elif np.abs(dy) < tol:
        x = np.linspace(x1, x2, nx)
        y = np.ones_like(x)*y1
    else:
        x = np.linspace(x1, x2, nx)
        slope = (y2-y1)/(x2-x1)
        y = slope*(x-x1)+y1
    return x, y


def plotMagSurvey2D(prob2D, susc, Einc, Edec, Bigrf, comp, irt,  Q, rinc, rdec, fig=None, axs1=None, axs2=None):

    import matplotlib.gridspec as gridspec

    # The MAG problem created is stored in result[1]
    # prob2D = Box.result[1]

    if fig is None:
        fig = plt.figure(figsize=(18*1.5,3.4*1.5))

        plt.rcParams.update({'font.size': 14})
        gs1 = gridspec.GridSpec(2, 7)
        gs1.update(left=0.05, right=0.48, wspace=0.05)

    if axs1 is None:
        axs1 = plt.subplot(gs1[:2, :3])

    if axs2 is None:
        axs2 = plt.subplot(gs1[0, 4:])

    axs1.axis("equal")

    prob2D.Bdec, prob2D.Binc, prob2D.Bigrf = Edec, Einc, Bigrf
    prob2D.Q, prob2D.rinc, prob2D.rdec = Q, rinc, rdec
    prob2D.uType, prob2D.mType = comp, 'total'
    prob2D.susc = susc

    # Compute fields from prism
    b_ind, b_rem = prob2D.fields()

    if irt == 'total':
        out = b_ind + b_rem

    elif irt == 'induced':
        out = b_ind

    else:
        out = b_rem

    X, Y = np.meshgrid(prob2D.survey.xr, prob2D.survey.yr)

    dat = axs1.contourf(X,Y, np.reshape(out, (X.shape)).T, 25)
    cb = plt.colorbar(dat, ax=axs1, ticks=np.linspace(out.min(), out.max(), 5))
    cb.set_label("nT")

    axs1.plot(X, Y, '.k')

    # Compute fields on the line by creating a similar mag problem
    prob2D.survey.profile
    if prob2D.survey.profile == "EW":
        x1, x2, y1, y2 = prob2D.survey.xr[0], prob2D.survey.xr[-1], 0., 0.
    elif prob2D.survey.profile == "NS":
        x1, x2, y1, y2 = 0., 0., prob2D.survey.yr[0], prob2D.survey.yr[-1]
    elif prob2D.survey.profile == "45N":
        x1, x2, y1, y2 = prob2D.survey.xr[0], prob2D.survey.xr[-1], prob2D.survey.yr[0], prob2D.survey.yr[-1]

    x, y = linefun(x1, x2, y1, y2, prob2D.survey.npts2D)
    xyz_line = np.c_[x, y, np.ones_like(x)*prob2D.survey.rx_h]
    # Create problem
    prob1D = MAG.problem()
    srvy1D = MAG.survey()
    srvy1D._rxLoc = xyz_line

    prob1D.prism = prob2D.prism
    prob1D.survey = srvy1D

    prob1D.Bdec, prob1D.Binc, prob1D.Bigrf = Edec, Einc, Bigrf
    prob1D.Q, prob1D.rinc, prob1D.rdec = Q, rinc, rdec
    prob1D.uType, prob1D.mType = comp, 'total'
    prob1D.susc = susc

    # Compute fields from prism
    out_linei, out_liner = prob1D.fields()

    #out_linei, out_liner = getField(p, xyz_line, comp, 'total')
    #out_linei = getField(p, xyz_line, comp,'induced')
    #out_liner = getField(p, xyz_line, comp,'remanent')

    out_linet = out_linei+out_liner

    distance = np.sqrt((x-x1)**2.+(y-y1)**2.)


    axs1.plot(x,y, 'w.', ms=3)

    axs1.text(x[0], y[0], 'A', fontsize=16, color='w')
    axs1.text(x[-1], y[-1], 'B', fontsize=16,
             color='w', horizontalalignment='right')

    axs1.set_xlabel('Easting (X; m)')
    axs1.set_ylabel('Northing (Y; m)')
    axs1.set_xlim(X.min(), X.max())
    axs1.set_ylim(Y.min(), Y.max())
    axs1.set_title(irt+' '+comp)

    axs2.plot(distance, out_linei, 'b.-')
    axs2.plot(distance, out_liner, 'r.-')
    axs2.plot(distance, out_linet, 'k.-')
    axs2.set_xlim(distance.min(), distance.max())

    axs2.set_xlabel("Distance (m)")
    axs2.set_ylabel("Magnetic field (nT)")

    axs2.text(distance.min(), out_linei.max()*0.8, 'A', fontsize = 16)
    axs2.text(distance.max()*0.97, out_linei.max()*0.8, 'B', fontsize = 16)
    axs2.legend(("induced", "remanent", "total"), bbox_to_anchor=(0.5, -0.3))
    axs2.grid(True)
    plt.show()

    return True


def fitline(Box):

    def profiledata(data, Binc, Bdec, Bigrf, x0, depth, susc, Q, rinc, rdec, update):

        prob = Box.result[1]
        prob.prism.z0 = -depth

        return plotProfile(prob, x0, data, Binc, Bdec, Bigrf, susc, Q, rinc, rdec)

    Q = widgets.interactive(profiledata, data=widgets.ToggleButtons(options=['MonSt','WedTA','WedSt']),\
             Binc=widgets.FloatSlider(min=-90.,max=90,step=5,value=90,continuous_update=False),\
             Bdec=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False),\
             Bigrf=widgets.FloatSlider(min=54000.,max=55000,step=10,value=54500,continuous_update=False),\
             x0=widgets.FloatSlider(min=5., max=25., step=0.1, value=15.), \
             depth=widgets.FloatSlider(min=0.,max=2.,step=0.05,value=0.5), \
             susc=widgets.FloatSlider(min=0., max=800.,step=5., value=1.),\
             Q=widgets.FloatSlider(min=0., max=10.,step=0.1, value=0.),\
             rinc=widgets.FloatSlider(min=-180., max=180.,step=1., value=0.),\
             rdec=widgets.FloatSlider(min=-180., max=180.,step=1., value=0.),
             update=widgets.ToggleButton(description='Refresh', value=False) \
             )
    return Q


def ViewMagSurvey2DInd(Box):


    def MagSurvey2DInd(susc, Einc, Edec, Bigrf, comp, irt, Q, rinc, rdec, update):

        # Get the line extent from the 2D survey for now
        prob = Box.result[1]

        return plogMagSurvey2D(prob, susc, Einc, Edec, Bigrf, comp, irt, Q, rinc, rdec)

    out = widgets.interactive (MagSurvey2DInd
                    ,susc=widgets.FloatSlider(min=0,max=200,step=0.1,value=0.1,continuous_update=False) \
                    ,Einc=widgets.FloatSlider(min=-90.,max=90,step=5,value=90,continuous_update=False) \
                    ,Edec=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False) \
                    ,Bigrf=widgets.FloatSlider(min=53000.,max=55000,step=10,value=54500,continuous_update=False) \
                    ,comp=widgets.ToggleButtons(options=['tf','bx','by','bz'])
                    ,irt=widgets.ToggleButtons(options=['induced','remanent', 'total'])
                    ,Q=widgets.FloatSlider(min=0.,max=10,step=1,value=0,continuous_update=False) \
                    ,rinc=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False) \
                    ,rdec=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False) \
                    ,update=widgets.ToggleButton(description='Refresh', value=False) \
                    )
    return out

def ViewPrismAndData():


    Q = widgets.interactive(PrismAndData \
                            , dx=widgets.FloatSlider(min=1e-4, max=500., step=10.0, value=120, continuous_update=False) \
                            , dy=widgets.FloatSlider(min=1e-4, max=500., step=10.0, value=120, continuous_update=False) \
                            , dz=widgets.FloatSlider(min=1e-4, max=100., step=1.0, value=7, continuous_update=False) \
                            , depth=widgets.FloatSlider(min=0., max=100., step=0.1, value=10, continuous_update=False)\
                            , pinc=(-90., 90., 5.) \
                            , pdec=(-90., 90., 5.) \
                            , susc=widgets.FloatSlider(min=0,max=10.0, value=5.0,step=0.5)
                            , npts2D=widgets.IntSlider(min=5, max=100, step=5, value=20, continuous_update=False) \
                            , xylim=widgets.FloatSlider(min=1, max=200, step=10, value=100, continuous_update=False) \
                            , rx_h=widgets.FloatSlider(min=0.1, max=2.5, step=0.1, value=1.5, continuous_update=False) \
                            , profile=widgets.ToggleButtons(options=['EW', 'NS', '45N'])
                            , View_elev=widgets.FloatSlider(min=-90, max=90, step=5, value=20, continuous_update=False) \
                            , View_azim=widgets.FloatSlider(min=0, max=360, step=5, value=250, continuous_update=False) \
                            , Binc=widgets.FloatSlider(min=-90.,max=90,step=5,value=45,continuous_update=False) \
                            , Bdec=widgets.FloatSlider(min=-90.,max=90,step=5,value=20,continuous_update=False) \
                            , Bigrf=widgets.FloatSlider(min=45000.,max=70000,step=100,value=50100,continuous_update=False)
                            )

    return Q


def PipeByPrism(depth, pdec, npts2D, xylim, rx_h, profile, View_elev, View_azim):
    #p = definePrism(dx, dy, dz, depth,pinc=pinc, pdec=pdec, susc = 1., Einc=90., Edec=0., Bigrf=1e-6)
    p = definePrism()
    p.dx=0.1
    p.dy=100
    p.dz=0.1
    p.y0=0.0
    p.z0 =-depth
    p.x0=0.0

    p.pinc, p.pdec = 0.0, pdec

    srvy = MAG.survey()
    srvy.rx_h, srvy.npts2D, srvy.xylim, srvy.profile = rx_h, npts2D, xylim, profile

    # Create problem
    prob = MAG.problem()
    prob.prism = p

#    p.x0=x0
    p.dec=pdec
    p.z0=-depth

    
    prob.survey = srvy

    return plotObj3D(p, rx_h, View_elev, View_azim, npts2D, xylim, profile=profile), prob

def ViewPipeByPrism(depth, dec, xylim=3.):
    elev, azim = 20, 250
    npts2D = 20
    Q=widgets.interactive(PipeByPrism \
                            , depth=widgets.FloatSlider(min=0., max=10., step=0.1, value=depth, continuous_update=False,description='Pipe z')\
                            , pdec=widgets.FloatSlider(min=-90., max=90., step=5.,value=90.0,description='Pipe dec') \
                            , npts2D=widgets.IntSlider(min=5, max=100, step=5, value=npts2D, continuous_update=False,description='npts grid') \
                            , xylim=widgets.FloatSlider(min=1, max=50, step=5, value=10.0, continuous_update=False,description='xy max') \
                            , rx_h=widgets.FloatSlider(min=0.1, max=2.5, step=0.1, value=2.0, continuous_update=False,description='survey height') \
                            , profile=widgets.ToggleButtons(options=['NS','EW',  '45N'])
                            , View_elev=widgets.FloatSlider(min=-90, max=90, step=5, value=0.0, continuous_update=False,description='view elev') \
                            , View_azim=widgets.FloatSlider(min=0, max=360, step=5, value=90.0, continuous_update=False,description='view azim')
                            )
    return Q


def PlotPipeField(prob2D,Binc,Bdec,Bigrf,susc,fig=None,axs=None, dolegend=True):


    
    Q=0.0
    rinc=0.0
    rdec=0.0

    p=prob2D.prism
    
    dx = 0.0#x0 - prob2D.survey.xylim

    prob2D.survey.profile
    if prob2D.survey.profile == "EW":
        x1, x2, y1, y2 = prob2D.survey.xr[0]-dx, prob2D.survey.xr[-1]-dx, 0., 0.
    elif prob2D.survey.profile == "NS":
        x1, x2, y1, y2 = 0., 0., prob2D.survey.yr[0]-dx, prob2D.survey.yr[-1]-dx
    elif prob2D.survey.profile == "45N":
        x1, x2, y1, y2 = prob2D.survey.xr[0], prob2D.survey.xr[-1], prob2D.survey.yr[0], prob2D.survey.yr[-1]

    x, y = linefun(x1, x2, y1, y2, 100)
    xyz_line = np.c_[x, y, np.ones_like(x)*prob2D.survey.rx_h]

    x0=np.sqrt((x1-p.x0)**2.+(y1-p.y0)**2.)


    distance = np.sqrt((x-x1)**2.+(y-y1)**2.)

    xlim = [distance[0],distance[-1]]
    prob1D = MAG.problem()
    srvy1D = MAG.survey()
    srvy1D._rxLoc = xyz_line

    prob1D.prism = p
    prob1D.survey = srvy1D

    prob1D.Bdec, prob1D.Binc, prob1D.Bigrf = Bdec, Binc, Bigrf
    prob1D.Q, prob1D.rinc, prob1D.rdec = Q, rinc, rdec
    prob1D.uType, prob1D.mType = 'tf', 'total'
    prob1D.susc = susc

    # Compute fields from prism
    magi, magr = prob1D.fields()

    print(np.max(magi),np.max(magr))

    if fig is None:
        f = plt.figure(figsize = (10, 5))
        gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])

    if axs is None:
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
    else:
        ax0=axs[0]
        ax1=axs[1]

#    x0plot=x0+prob2D.survey.xr[0]
        




    ax0.set_ylabel("Total field anomaly (nT)")

    ax0.grid(True)

    if prob2D.survey.profile == "EW":

        if ax1 is not 0:
            ax1.set_xlabel("Easting (m)")
            ax1.plot(x0+prob2D.survey.xr[0], p.z0, 'ko',label='pipe')
            ax1.plot(xlim+prob2D.survey.xr[0], [prob2D.survey.rx_h, prob2D.survey.rx_h], 'b--',label='Magnetometer Height')
        if ax0 is not 0:
            ax0.plot(distance+prob2D.survey.xr[0], magi+magr, 'k', label='total')
            ax0.set_xlabel("Easting (m)")

    elif prob2D.survey.profile == "NS":
        if ax1 is not 0:
            ax1.set_xlabel("Northing (m)")
            ax1.plot(x0+prob2D.survey.yr[0], p.z0, 'ko',label='pipe')
            ax1.plot(xlim+prob2D.survey.yr[0], [prob2D.survey.rx_h, prob2D.survey.rx_h], 'b--',label='Magnetometer Height')
        if ax0 is not 0:
            ax0.plot(distance+prob2D.survey.yr[0], magi+magr, 'k', label='total')
            ax0.set_xlabel("Northing (m)")

    elif prob2D.survey.profile == "45N":
        if ax1 is not 0:
            ax1.set_xlabel("Distance SW-NE (m)")
            ax1.plot(x0, p.z0, 'ko',label='pipe')
            ax1.plot(xlim, [prob2D.survey.rx_h, prob2D.survey.rx_h], 'b--',label='Magnetometer Height')
        if ax0 is not 0:
            ax0.plot(distance, magi+magr, 'k', label='total')
            ax0.set_xlabel("Distance SW-NE (m)")

    if ax1 is not 0:
        ax1.grid(True)
        if dolegend is True:
            ax1.legend(fontsize=10)
        ax1.set_ylim([-10,5])

            # ax1.invert_yaxis()

#    ax1.text(x0+0.5, p.z0, 'Pipe', color='k')
#    ax1.text(xlim[0]+1.,prob2D.survey.rx_h-0.5, 'Magnetometer height', color='b')


    plt.tight_layout()
    plt.show()

    return True
def PlotPipeField_withQ(prob2D,Binc,Bdec,Bigrf,susc,Q,rinc,rdec,fig=None,axs=None, dolegend=True):


    p=prob2D.prism
    
    dx = 0.0#x0 - prob2D.survey.xylim

    prob2D.survey.profile
    if prob2D.survey.profile == "EW":
        x1, x2, y1, y2 = prob2D.survey.xr[0]-dx, prob2D.survey.xr[-1]-dx, 0., 0.
    elif prob2D.survey.profile == "NS":
        x1, x2, y1, y2 = 0., 0., prob2D.survey.yr[0]-dx, prob2D.survey.yr[-1]-dx
    elif prob2D.survey.profile == "45N":
        x1, x2, y1, y2 = prob2D.survey.xr[0], prob2D.survey.xr[-1], prob2D.survey.yr[0], prob2D.survey.yr[-1]

    x, y = linefun(x1, x2, y1, y2, 100)
    xyz_line = np.c_[x, y, np.ones_like(x)*prob2D.survey.rx_h]

    x0=np.sqrt((x1-p.x0)**2.+(y1-p.y0)**2.)


    distance = np.sqrt((x-x1)**2.+(y-y1)**2.)

    xlim = [distance[0],distance[-1]]
    prob1D = MAG.problem()
    srvy1D = MAG.survey()
    srvy1D._rxLoc = xyz_line

    prob1D.prism = p
    prob1D.survey = srvy1D

    prob1D.Bdec, prob1D.Binc, prob1D.Bigrf = Bdec, Binc, Bigrf
    prob1D.Q, prob1D.rinc, prob1D.rdec = Q, rinc, rdec
    prob1D.uType, prob1D.mType = 'tf', 'total'
    prob1D.susc = susc

    # Compute fields from prism
    magi, magr = prob1D.fields()

    if fig is None:
        f = plt.figure(figsize = (10, 5))
        gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])

    if axs is None:
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
    else:
        ax0=axs[0]
        ax1=axs[1]

    if ax1 is not 0:
        ax1.plot(x0, p.z0, 'ko',label='pipe')
        ax1.plot(xlim, [prob2D.survey.rx_h, prob2D.survey.rx_h], 'b--',label='Magnetometer Height')
        ax1.set_xlabel("Northing (m)")
        ax1.set_ylabel("Depth (m)")
        ax1.grid(True)
        if dolegend is True:
            ax1.legend(fontsize=10)
        ax1.set_ylim([-10,5])


    ax0.plot(distance, magi+magr, 'k', label='total')


    ax0.set_ylabel("Total field anomaly (nT)")

    ax0.grid(True)

    if prob2D.survey.profile == "EW":
        if ax1 is not 0:
            ax1.set_xlabel("Easting (m)")
        ax0.set_xlabel("Easting (m)")

    elif prob2D.survey.profile == "NS":
        if ax1 is not 0:
            ax1.set_xlabel("Northing (m)")
        ax0.set_xlabel("Northing (m)")

    elif prob2D.survey.profile == "45N":
        if ax1 is not 0:
            ax1.set_xlabel("Distance SW-NE (m)")
        ax0.set_xlabel("Distance SW-NE (m)")

            # ax1.invert_yaxis()

#    ax1.text(x0+0.5, p.z0, 'Pipe', color='k')
#    ax1.text(xlim[0]+1.,prob2D.survey.rx_h-0.5, 'Magnetometer height', color='b')


    plt.tight_layout()
    plt.show()

    return True



def ViewPipeIncField(Box):

    prob=Box.result[1]

    def SetupPipeField(Binc,Bdec,BF,susc,rx_h,depth):

        prob.survey.rx_h=rx_h
        prob.prism.z0=-depth
        prob.prism.x0=0.0
#        prob.survey.profile=profile

        return PlotPipeField(prob,Binc,Bdec,BF,susc)



    Q = widgets.interactive(SetupPipeField \
                          , Binc=widgets.FloatSlider(min=-90.,max=90,step=5,value=65,continuous_update=False) \
                          , Bdec=widgets.FloatSlider(min=-90.,max=90,step=5,value=-20,continuous_update=False) \
                          , BF=widgets.FloatSlider(min=45000.,max=70000,step=100,value=50100,continuous_update=False) \
                          , susc=widgets.FloatSlider(min=0.,max=200,step=0.1,value=100,continuous_update=False) \
                          , rx_h=widgets.FloatSlider(min=0.1, max=2.5, step=0.1, value=prob.survey.rx_h, continuous_update=False,description='survey height') \
                          , depth=widgets.FloatSlider(min=0., max=10., step=0.1, value=-prob.prism.z0, continuous_update=False,description='Pipe z') \
#                          , profile=prob.survey.profile
                          )

    return Q

def ViewPipeIncField_withQ(Box):

    prob=Box.result[1]

    def SetupPipeField_withQ(Binc,Bdec,BF,susc,Qr,rinc,rdec,rx_h,depth,x0,profile):
        prob.survey.rx_h=rx_h
        prob.prism.z0=-depth
        prob.prism.x0=x0
        prob.survey.profile=profile

        return PlotPipeField_withQ(prob,Binc,Bdec,BF,susc,Qr,rinc,rdec)


    Q = widgets.interactive(SetupPipeField_withQ \
                          , Binc=widgets.FloatSlider(min=-90.,max=90,step=5,value=45,continuous_update=False) \
                          , Bdec=widgets.FloatSlider(min=-90.,max=90,step=5,value=20,continuous_update=False) \
                          , BF=widgets.FloatSlider(min=45000.,max=70000,step=100,value=50100,continuous_update=False) \
                          , susc=widgets.FloatSlider(min=0.,max=200,step=0.1,value=100,continuous_update=False)
                          , Qr=widgets.FloatSlider(min=0.,max=1.0,step=0.1,value=0.0,continuous_update=False)
                          , rinc=widgets.FloatSlider(min=-90.,max=90.,step=1.0,value=0.0,continuous_update=False)
                          , rdec=widgets.FloatSlider(min=-90.,max=90.,step=1.0,value=100,continuous_update=False)
                          , rx_h=widgets.FloatSlider(min=0.1, max=2.5, step=0.1, value=prob.survey.rx_h, continuous_update=False) \
                          , depth=widgets.FloatSlider(min=0., max=10., step=0.1, value=prob.prism.z0, continuous_update=False)
                          , x0=widgets.FloatSlider(min=-5.0, max=5.0, step=0.1, value=prob.prism.z0, continuous_update=False)
                          , profile=widgets.ToggleButtons(options=['EW', 'NS', '45N']) \
                          )

    return Q

def PlotNPrisms(x,y,z,dx,dy,dz,pdec,pinc, nprisms, npts2D, xylim, rx_h, profile, View_elev, View_azim):
    #p = definePrism(dx, dy, dz, depth,pinc=pinc, pdec=pdec, susc = 1., Einc=90., Edec=0., Bigrf=1e-6)

    fig=plt.figure()
    axs=fig.add_subplot(111,projection='3d')
    
    cstring=['c','w','m','y','g','r','b']

    problist=[]


    for ii in range(nprisms):
        p = definePrism()
        p.dx=dx[ii]
        p.dy=dy[ii]
        p.dz=dz[ii]
        p.x0=x[ii]
        p.y0=y[ii]
        p.z0 =-z[ii]
        p.pinc=pinc[ii]
        p.pdec =pdec[ii]
        
        srvy = MAG.survey()
        srvy.rx_h, srvy.npts2D, srvy.xylim, srvy.profile = rx_h, npts2D, xylim, profile


        # Create problem
        prob = MAG.problem()
        prob.prism = p

        prob.survey = srvy

        problist.append(prob)

        plotObj3D(p, rx_h, View_elev, View_azim, npts2D, xylim, profile=profile,fig=fig,axs=axs,cstring=cstring[ii])

    plt.show()
    return problist

def View3Prisms(npts2D=20,xylim=100,rx_h=1.9,profile='NS',elev=20,azim=250):

    def Plot3Prisms(x1,x2,x3,y1,y2,y3,z1,z2,z3,dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,pd1,pd2,pd3,pi1,pi2,pi3, npts2D, xylim, rx_h, profile, View_elev, View_azim):
        nprisms=3
        x=np.array([x1,x2,x3])
        y=np.array([y1,y2,y3])
        z=np.array([z1,z2,z3])
        dx=np.array([dx1,dx2,dx3])
        dy=np.array([dy1,dy2,dy3])
        dz=np.array([dz1,dz2,dz3])
        pdec=np.array([pd1,pd2,pd3])
        pinc=np.array([pi1,pi2,pi3])

        return PlotNPrisms(x,y,z,dx,dy,dz,pdec,pinc,nprisms,npts2D,xylim,rx_h,profile,View_elev,View_azim)

    Q = widgets.interactive(Plot3Prisms \
                            , x1=widgets.FloatSlider(min=-100., max=100., step=1, value=1.0, continuous_update=False)\
                            , y1=0.0 \
                            , z1=widgets.FloatSlider(min=0., max=100., step=0.1, value=1.0, continuous_update=False)\
                            , dx1=widgets.FloatSlider(min=0., max=500., step=0.1, value=1.0, continuous_update=False)\
                            , dy1=widgets.FloatSlider(min=0., max=500., step=0.1, value=1.0, continuous_update=False)\
                            , dz1=widgets.FloatSlider(min=0., max=100., step=0.1, value=1.0, continuous_update=False)\
                            , pd1=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , pi1=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , x2=widgets.FloatSlider(min=-10., max=10., step=0.1, value=1.0, continuous_update=False)\
                            , y2=0.0 \
                            , z2=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dx2=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dy2=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dz2=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , pd2=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , pi2=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , x3=widgets.FloatSlider(min=-10., max=10., step=0.1, value=1.0, continuous_update=False)\
                            , y3=0.0 \
                            , z3=widgets.FloatSlider(min=0., max=10., step=0.1, value=1.0, continuous_update=False)\
                            , dx3=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dy3=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dz3=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , pd3=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , pi3=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , npts2D=widgets.IntSlider(min=5, max=100, step=5, value=npts2D, continuous_update=False) \
                            , xylim=widgets.FloatSlider(min=1, max=100, step=10, value=xylim, continuous_update=False) \
                            , rx_h=widgets.FloatSlider(min=0.1, max=2.5, step=0.1, value=rx_h, continuous_update=False) \
                            , profile=widgets.ToggleButtons(options=['EW', 'NS', '45N'])
                            , View_elev=widgets.FloatSlider(min=-90, max=90, step=5, value=elev, continuous_update=False) \
                            , View_azim=widgets.FloatSlider(min=0, max=360, step=5, value=azim, continuous_update=False)
                            )

    return Q



def Gpoly(xobs,zobs,xcorn,zcorn,rho):
#From Blakely's book, computes the vertical component of gravity for a 2D body with a polygonal cross-section. 
#"Axes are right-handed with y-axis parallel to long direction of body and z axis vertical down, observation 
#points are xobs, zobs; xcorn, zcorn are the corners of the body arranged in clockwise order when viewed with 
#x axis to the right.  Density is rho in kg/m^3" output is in mgals?


    ncorn=np.size(xcorn)

    if ncorn is not np.size(zcorn):
        print("error, number of xcorn and zcorn points must be the same")

    if np.min(ncorn)<0.0:
        print("one of your corners is above the ground!")

    gamma=6.67e-11
    si2mg=1.0e5
    km2m=1000.0
    thresh=1.0e-8

    g=np.zeros_like(xobs)
    
    for ii in range(ncorn):
        if ii is ncorn-1:
            i2=0
        else:
            i2=ii+1
        
        x1=xcorn[ii]-xobs
        z1=zcorn[ii]-zobs

        x2=xcorn[i2]-xobs
        z2=zcorn[i2]-zobs

        r1sq=x1**2+z1**2
        r2sq=x2**2+z2**2

        if np.min(r1sq)<thresh or np.min(r2sq)<thresh:
            print("error one polynomial point is at the observation point")

        denom=z2-z1

        if denom.any()<thresh:
            for ix in range(np.size(xobs)):
                if denom[ix]<thresh:
                    denom[ix]=thresh

        alpha=(x2-x1)/denom

        beta=(x1*z2-x2*z1)/denom

        factor=beta/(1+alpha**2.)
        
        term1=0.5*(np.log(r2sq)-np.log(r1sq))
        term2=np.arctan2(z2,x2)-np.arctan2(z1,x1)

        g=g+factor*(term1-alpha*term2)

    g=g*2.0*rho*gamma*si2mg*km2m

    return g


def ComputeGravData(problist,rhovec,x,z):
    
    grav=np.zeros_like(x)

    plt.figure()

    for prob,rho in zip(problist,rhovec):
        p=prob.prism
        
        xvec,zvec=Get2DPrisms(p)

        print(xvec,zvec)

        noise=np.random.normal(0.1,2.0,np.size(grav))
        grav=grav+Gpoly(x,z,xvec,-zvec,rho)+noise

        plt.plot(xvec,zvec)
    
    plt.show()    
    return grav



def ExtractLog(xvec,zvec,x,z,val1,val2):
    
    vals=val1*np.ones_like(z)

    if x>np.max(xvec) or x < np.min(xvec):
        return vals
    
    slope=(zvec[1]-zvec[0])/(xvec[1]-xvec[0])
    int1=zvec[0]+slope*(x-xvec[0])
    int2=zvec[2]+slope*(x-xvec[2])
    
    i1=np.argmin(np.abs(z-int1))
    i2=np.argmin(np.abs(z-int2))
    
    vals[i1:i2]=val2
    
    return vals
    

def Get2DPrisms(p):

    x1, x2 = p.xn[0]-p.xc, p.xn[1]-p.xc
    y1, y2 = p.yn[0]-p.yc, p.yn[1]-p.yc
    z1, z2 = p.zn[0]-p.zc, p.zn[1]-p.zc
    block_xyz = np.asarray([[x1, x1, x2, x2, x1, x1, x2, x2],
                           [y1, y2, y2, y1, y1, y2, y2, y1],
                           [z1, z1, z1, z1, z2, z2, z2, z2]])

    R = Utils.rotationMatrix(p.pinc, p.pdec)

    xyz = R.dot(block_xyz).T

    x1r=xyz[0,0]
    x2r=xyz[1,0]
    x3r=xyz[4,0]
    x4r=xyz[5,0]
    z1r=xyz[0,2]
    z2r=xyz[1,2]
    z3r=xyz[4,2]
    z4r=xyz[5,2]

    xvec=[x3r,x4r,x2r,x1r,x3r]+p.xc
    zvec=[z3r,z4r,z2r,z1r,z3r]+p.zc

    return xvec,zvec



def ViewPrism(npts2D=20,xylim=100,rx_h=1.9,profile='NS',elev=20,azim=250):

    Q = widgets.interactive(Prism \
                            , x1=widgets.FloatSlider(min=-100., max=100., step=1, value=1.0, continuous_update=False)\
                            , y1=0
                            , z1=widgets.FloatSlider(min=0., max=100., step=0.1, value=1.0, continuous_update=False)\
                            , dx1=widgets.FloatSlider(min=0., max=500., step=0.1, value=1.0, continuous_update=False)\
                            , dy1=widgets.FloatSlider(min=0., max=500., step=0.1, value=1.0, continuous_update=False)\
                            , dz1=widgets.FloatSlider(min=0., max=100., step=0.1, value=1.0, continuous_update=False)\
                            , pd1=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , pi1=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , x2=widgets.FloatSlider(min=-10., max=10., step=0.1, value=1.0, continuous_update=False)\
                            , y2=0.0 \
                            , z2=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dx2=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dy2=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dz2=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , pd2=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , pi2=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , x3=widgets.FloatSlider(min=-10., max=10., step=0.1, value=1.0, continuous_update=False)\
                            , y3=0.0 \
                            , z3=widgets.FloatSlider(min=0., max=10., step=0.1, value=1.0, continuous_update=False)\
                            , dx3=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dy3=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , dz3=widgets.FloatSlider(min=0., max=20., step=0.1, value=1.0, continuous_update=False)\
                            , pd3=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , pi3=widgets.FloatSlider(min=-90., max=90., step=5.,value=0.0) \
                            , npts2D=widgets.IntSlider(min=5, max=100, step=5, value=npts2D, continuous_update=False) \
                            , xylim=widgets.FloatSlider(min=1, max=100, step=10, value=xylim, continuous_update=False) \
                            , rx_h=widgets.FloatSlider(min=0.1, max=2.5, step=0.1, value=rx_h, continuous_update=False) \
                            , profile=widgets.ToggleButtons(options=['EW', 'NS', '45N'])
                            , View_elev=widgets.FloatSlider(min=-90, max=90, step=5, value=elev, continuous_update=False) \
                            , View_azim=widgets.FloatSlider(min=0, max=360, step=5, value=azim, continuous_update=False)
                            )

    return Q



def FitGravMagData():

    def PlotGravMagData(width1,thick1,xpos1,depth1,susc1,rho1,width2,thick2,xpos2,depth2,susc2,rho2,width3,thick3,xpos3,depth3,susc3,rho3,pinc):
        
     
        data=np.loadtxt('Lab10_Data.txt', encoding = "utf-16")
        x=data[:,0]
        true_grav=data[:,2]
        true_mag=data[:,3]

        
        dx=500.0
        pdec=90.0
        problist=[]
        xylim=100.
        npts2D=20
        rx_h=1.9
        Binc=79.
        Bdec=-11.2
        Bigrf=57076.
        profile='EW'
        
        dy=np.array([width1,width2,width3])
        dz=np.array([thick1,thick2,thick3])
        x=np.array([xpos1,xpos2,xpos3])
        y=0.0
        z=-np.array([depth1,depth2,depth3])
       
        suscvec=np.array([susc1,susc2,susc3])
        rhovec=np.array([rho1,rho2,rho3])
        
        Q=0.0
        rinc=0.0
        rdec=0.0
    
        linex=np.linspace(0,1000,266)
        liney=np.zeros_like(linex)
        linez=np.ones_like(linex)*rx_h
        linez_grav=np.zeros_like(linex)
        
        line_xyz=np.c_[linex,liney,linez]
        
        Magdat=np.zeros_like(linex)
        Gravdat=np.zeros_like(linex)
        
        fig,ax=plt.subplots(3, sharex=True)
        fig.subplots_adjust(hspace=0)

#        [:-1]]
        plt.rcParams.update({'font.size': 10})
        
        ax_grav=ax[2]
        ax_mag=ax[1]
        ax_prism=ax[0]
        
        ax_grav=plt.subplot(3,1,1)
        ax_mag=plt.subplot(3,1,2)
        ax_prism=plt.subplot(3,1,3)
        
        ax_grav.set_ylabel('Gravity Anomaly (mGal)',fontsize=10)

        ax_mag.set_ylabel('Total Field (nT)')

        ax_prism.set_ylabel('Depth (m)',fontsize=10)
        ax_prism.set_xlabel('Easting (m)',fontsize=10)
        
        ax_grav.plot(linex,true_grav)
        ax_mag.plot(linex,true_mag,label='data')
        
        for ii in range(3):
            #setup problem
            prob=MAG.problem()
            p=definePrism()
            p.dx=dx
            p.dy=dy[ii]
            p.dz=dz[ii]
            p.x0=x[ii]
            p.y0=0.0
            p.z0=z[ii]
            p.pinc=pinc
            p.pdec=90.0
            
            srvy=MAG.survey()
            srvy.rx_h, srvy.npts2D, srvy.xylim, srvy.profile = rx_h, npts2D, xylim, profile
            
            srvy._rxLoc=line_xyz
            
            prob.prism=p
            prob.survey=srvy
            
            prob.Bdec,prob.Binc,prob.Bigrf=Bdec,Binc,Bigrf
            
            prob.susc=suscvec[ii]
            
            prob.Q, prob.rinc, prob.rdec = Q, rinc, rdec
            prob.uType, prob.mType = 'tf', 'total'
            
            #calculate data
            
            dati, datr=prob.fields()
            
            Magdat=Magdat+dati
            
            xvec,zvec=Get2DPrisms(p)
            
            Gravdat=Gravdat+Gpoly(linex,linez_grav,xvec,-zvec,rhovec[ii])
            
            ax_prism.plot(xvec,zvec)

        #do plotting
            
        ax_grav.plot(linex,Gravdat)
        ax_mag.plot(linex,Magdat,label='modelled')
        ax_mag.yaxis.tick_right()
        ax_mag.yaxis.set_label_position("right")


    
        ax_prism.set_xlim([linex[0],linex[-1]])
        ax_grav.set_xlim([linex[0],linex[-1]])
        ax_mag.set_xlim([linex[0],linex[-1]])
        ax_prism.set_ylim([-60,0])
        
        ax_mag.legend(fontsize=10)
        
        ax_mag.grid(which="both")
        ax_grav.grid(which="both")
        ax_prism.grid(which="both")


        plt.setp([ax_grav.get_xticklabels()], visible=False)
        plt.setp([ax_mag.get_xticklabels()], visible=False)

        
        plt.show()
        
        return True
    Q = widgets.interactive(PlotGravMagData \
                            , xpos1=widgets.FloatSlider(min=-100, max=1000., step=10.0, value=500, continuous_update=False) \
                            , xpos2=widgets.FloatSlider(min=-100, max=1000., step=10.0, value=500, continuous_update=False) \
                            , xpos3=widgets.FloatSlider(min=-100, max=1000., step=10.0, value=500, continuous_update=False) \
                            , width1=widgets.FloatSlider(min=1e-4, max=750., step=10.0, value=20, continuous_update=False) \
                            , thick1=widgets.FloatSlider(min=1e-4, max=100., step=1.0, value=20, continuous_update=False) \
                            , depth1=widgets.FloatSlider(min=0., max=100., step=0.1, value=-10, continuous_update=False)\
                            , width2=widgets.FloatSlider(min=1e-4, max=750., step=10.0, value=20, continuous_update=False) \
                            , thick2=widgets.FloatSlider(min=1e-4, max=100., step=1.0, value=20, continuous_update=False) \
                            , depth2=widgets.FloatSlider(min=0., max=100., step=0.1, value=-10, continuous_update=False)\
                            , width3=widgets.FloatSlider(min=1e-4, max=750., step=10.0, value=20, continuous_update=False) \
                            , thick3=widgets.FloatSlider(min=1e-4, max=100., step=1.0, value=20, continuous_update=False) \
                            , depth3=widgets.FloatSlider(min=0., max=100., step=0.1, value=-10, continuous_update=False)\
                            , pinc=(-90., 90., 5.) \
                            , susc1=widgets.FloatSlider(min=0,max=1.0, step=0.01)
                            , susc2=widgets.FloatSlider(min=0,max=1.0, step=0.01)
                            , susc3=widgets.FloatSlider(min=0,max=1.0, step=0.01)
                            , rho1=widgets.FloatSlider(min=-1,max=1.0, step=0.05)
                            , rho2=widgets.FloatSlider(min=-1,max=1.0, step=0.05)
                            , rho3=widgets.FloatSlider(min=-1,max=1.0, step=0.05)
                            )

    return Q

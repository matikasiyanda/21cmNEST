#!/usr/bin/env python

"""This is straightline.py
Jonathan Zwart
January 2016
"""
# import sys

#sys.path.insert(0,"/home/siyanda/siyanda/lib64/python2.6/site-packages/scipy/ ")
from mpi4py import MPI
import os,sys
import importlib
import numpy
import numpy as np
from math import pi,exp,log,sqrt,log
import pymultinest
#from scipy.special import erf
from priors import Priors
#import siyanda_functions as sm
#import siyanda1_fun as sm
import Absorption_S as AB
#import model_1 as sm
#import non_siyanda as sm
from scipy import integrate
#from profile_support import profile
from utils import medianArray

#--------------------------------------
import errno
def DIR(dirname):
 try:
  os.mkdir(dirname)
 except OSError as exc:
  if exc.errno != errno.EEXIST:
   raise exc
  pass
#------------------------------------



# Constants
sqrtTwo=sqrt(2.0)

if __name__ == '__main__' and len(sys.argv) < 5:
    print 'usage:'
    print
    print 'with MPI [NB mpiEXEC on some machines]:'
    print '       mpiexec -n NPROC ./straightline.py settings.py zg L0 ag IDnum  '
    print
    sys.exit(0)

param_file=sys.argv[-4]
try:
    execfile(param_file)
except IOError:
    from settings import *

#-------------------------------------------------------------------------------
#Inited=True

def chi(z):
    H0=68.
    oh_m=0.27
    oh_r=9.2e-5
    c=3.0e5
    a =np.sqrt(oh_r*(1.+z)**4. +oh_m*(1.+z)**3.+(1.-oh_m-oh_r))
    b= (H0)*a
    return (c/b)

def r(z):
    return integrate.quad(chi,0,z)[0]

def Rnter(f):
        """
        This an interplotion function for r(z) the comoving distnace
        units : Mpc
        bb: is an array of redshifts from 6 up to 20 with stepsize dz=0.1
        cc: contains the corresponding values of r(z) in bb
        both bb,cc have length 141
        """
        s=np.interp(f,bb,cc)
        return s


def Snu_Model_A(x,B,K,z_g,a,L0):
    """
    x: Observer's frequency [MHz]
    B: 1st fitting parameter for 21 Absorption feature
    K: 2st fitting parameter for 21 Absorption feature
    z_g: Emitting source's redshift position
    a: Source spectral index
    L0: Source's Luminosity

    observer frame
    xrest/xobs=1+z
    x=freq_obs
    z_g=z_source
    freq_break_
    """
    nu01 =x*1.0e-3
    var1=0.214*(L0/1.0)*(46.77)/( (1.0+z_g)**(1.+a) )
    var2 =(9142.7/ r(z_g) )**2.0
    var3 = var1*var2*(1.4/nu01)**a
    zo=(1420.4/x)-1.0

    if zo<z_g:#Absorption feature
        var=  B*np.log(1.0+zo) + K
        return  (1.0-np.exp(var) )*var3
    else:# no Absorption feature
        return 1.0*var3



def Int_S(x,x_i,B,K,z_g,a,L0):
    var=integrate.quad(lambda nu:Snu_Model_A(nu,B,K,z_g,a,L0),x,x_i)
    return var[0]/(x_i-x)


def Tsys(nu):
    '''
    Tsys: Represents the sky temperature at a given frequency
    T_sky: Sky temperature :units [k ]

    args:
    The only argument required is frequncy (nu)

    '''
    T_sky = 60.0*(300./nu)**(2.55)
    T=(0.1*T_sky+40.0)+( T_sky)
    return T



def sig_s(nu,dnu,t_int):
    var1= Tsys(nu)*2.0*1.38e-23
    if nu>110.:
        A_tt = 0.5*866.0*(925.0)*(110./nu)**2.0
    else:
        A_tt =0.5*866.0*(925.0)
    #A_tt = 0.5*866.0*(925.0)*(110./nu)**2.0
    var2= var1*1.0e26/ ( A_tt*np.sqrt(dnu*1.0e6*t_int*60.0*60.0) )
    return var2*1.0e3 #!!NB the 10 dividing is to improve the signal to noise



pri=Priors()
def myprior(cube, ndim, nparams):


    cube[0]=pri.GeneralPrior(cube[0],'U',5.0,20.0)#zfit
    cube[1]=pri.GeneralPrior(cube[1],'U',0.0,100.0)#a0=B
    cube[2]=pri.GeneralPrior(cube[2],'U',-100.0,00.0)#a1=K
    cube[3]=pri.GeneralPrior(cube[3],'LOG',1.0e-6,1.0e6)#L0
    cube[4]=pri.GeneralPrior(cube[4],'U',0.0,10.0)#sigma
    cube[5]=pri.GeneralPrior(cube[5],'U',0.0,2.0)#alpha



    return

#-------------------------------------------------------------------------------



notInited=True
def myloglike(cube, ndim, nparams):
    global yi,x,xm,y,notInited,sg,bb,cc
    if notInited:

        zin,Li,ai=float(sys.argv[-4]),float(sys.argv[-3]),float(sys.argv[-2])
        F=AB.FLUX(zin,Li,ai)
        bb,cc = zip(* AB.Memoize(r).memo.items())
        y=F.Result()[1]#
        x=F.Result()[0]#
        xm=medianArray(x)
        sg=F.Result()[2]#np.array( [ sig_s(i,dnu,T_int) for i in xm[:-1] ]  )
        np.random.seed(1234)
        yi=y+np.random.normal(0.0,sg)#input data + random.noramal noise
        notInited=False


    loglike=0.0

    for idatum in range(len(xm[:-1])):
        data=yi[idatum]
        zfit=cube[0]
        B=cube[1]
        K=cube[2]
        amp=cube[3]
        a=cube[5]
        model=Int_S(xm[idatum],xm[idatum+1],B,K,zfit,a,amp)
        if model < 0.0:
            print '+',
            return -1.0e99
        #print model
        sigma=cube[4]*sg[idatum]
        chisq = 0.5*((data-model)/sigma)**2.0
        prefactor=0.5*log(2.0*pi*sigma**2.0)
        loglike -= prefactor + chisq
    return loglike

#-------------------------------------------------------------------------------

def main():

    """
    """

    if not os.path.exists(outdir): DIR(outdir)

    n_params=6
    pymultinest.run(myloglike,myprior,n_params,resume=False,verbose=True,\
                multimodal=multimodal,max_modes=max_modes,write_output=True,\
                n_live_points=n_live_points,\
                evidence_tolerance=evidence_tolerance,\
                mode_tolerance=-1e90,seed=seed,max_iter=max_iter,\
                importance_nested_sampling=do_INS,\
                outputfiles_basename=os.path.join(outdir,outstem),init_MPI=False)

#contour_plot.contourTri(pylab.loadtxt('chains_test/1-post_equal_weights.dat'),line=True,outfile='chains_test/test.png',col=('red','blue'),labels=parameters,binsize=50,truth=plotTruth,ranges=plotRanges,autoscale=True)

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

#!/usr/bin/env python

"""
This is plot_tri.py
Jonathan Zwart
February 2016

Usage:

./plot_tri.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
from matplotlib import pyplot as plt
import pylab
#from profile_support import profile
from utils import *
import contour_plot

#param_file=sys.argv[-1]
#setf='%s.settings' % param_file
outdir=sys.argv[-5]

#-------------------------------------------------------------------------------

def main():

    """
    """

    #print 'Settings file is %s' % param_file

    # Import the settings variables
    #set_module=importlib.import_module(setf)
    #globals().update(set_module.__dict__)

    chain=pylab.loadtxt('%s/1-post_equal_weights.dat'%outdir)

    # Plot for publication
    line=True
    autoscale=True
    title=''
    extn='pdf'
    numbins=20
    nparams=chain.shape[-1]-1
    parameters=[p for p in range(nparams)]
    #labelDict={'0':r'm','1':r'c','2':r'sig'}
    #truth={'0':6.1,'1':-21.3,'2':1}
    #labelDict={0:r"$\alpha$",1:r"$L_{0}$",2:r"$ z_g$",3:r"$ \sigma$"}
    labelDict={0:r" $ z_g$",1:r"$a_{1}$",2:r"$a_{0}$",3:r"$ \frac{ L_{0} }{10^{26}[\rm{W/Hz} ] }  $",4:r" $  \sigma_{0}   $ ",5:r"$ \alpha $"}
    truth={0:14.0  ,1:12.41,2:-36.18,3:float(sys.argv[-3]) ,4:1.0,5:0.70}
    print float(sys.argv[-4]),sys.argv[-3]
    #ruth=None
    plotRanges={0:[10.0,18.0],1:[ 0.0,100.0],2:[-100.0,-10.0],3:[0.0,float(sys.argv[-2])*2.0], 4:[0.00,5.00],5:[0.5,0.9]}
    plotRanges=None
    #plotRanges={'0':[5.0,7.0],'1':[-20,-22.0],'2':[0.0,1.0]}
    bundle=contour_plot.contourTri(chain,\
                            line=line,\
                            outfile='%s/triangle.%s'%(outdir,extn),\
                            autoscale=autoscale,\
                            labels=parameters,truth=truth,labelDict=labelDict,\
                            binsize=numbins,\
                            ranges=plotRanges)
#                            ranges=plotRanges,truth=truth,\
#                            autoscale=autoscale,title=title,\
#                            binsize=binsize,labelDict=labelDict)

    #stats=fetchStats(outdir,parameters,plotTruth)
    #printLaTeX(parameters,stats,dump=outdir)


#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

"""
Parameter file for straightline.py
"""
#-------------------------------------------------------------------------------
import numpy as np
import os,errno,sys
outdir='chains_Mario_'+str(sys.argv[-1]) #+str(np.random.randint(0,1001))

#outdir=
#newpath='chains_Mario
#outdir=newpath#safe_make_folder(newpath)
#-------------------------------------------------------------------------------

# Set some MultiNEST parameters here
#n_live_points=1000
multimodal=False
max_modes=1
seed=1234 # [-1 for clock]

# Switch to INS
do_INS=False
n_live_points=500
#n_live_points=1000
outstem='1-'

max_iter=0
evidence_tolerance=0.36

#random.seed(1234)
	

#-------------------------------------------------------------------------------

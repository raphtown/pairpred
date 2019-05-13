# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 01:41:22 2012
@author: fayyaz
This module saves the myPDB objects for the pdb files saved in pdbpath in the 
present directory
These files are then used in dbdscrpp2.py
How to run: mpirun -np 8 python getPDBPKL.py
Note: if there is an error in processing a PDB file, the output is not saved
"""
from myPDB import *
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
import glob
import numpy as np


pdbpath='./DBD4/'
profilepath=pdbpath+'Profile/'
opath='./DBD4PKLN/'
f=glob.glob(pdbpath+'*.pdb')
N=len(f)*1.0
bindx=int(np.floor(rank*N/size))
eindx=int(np.floor((rank+1)*N/size))
for e in range(bindx,eindx):
    print str(e)
    fn=f[e]
    try:
        L=myPDB(fn,profilepath)
        L.save(bdir=opath)
    except:
        print "ERRRORRRRRR in "+fn+" # "+str(e)
        
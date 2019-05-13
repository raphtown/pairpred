# -*- coding: utf-8 -*-
"""
Created on DEC31 2012

@author: fayyaz
Descripton: Performs nested cross validation over a kernel file (*.dbK.pkl).
It chooses the optimal C value for each of the N folds fold by performing N-1
fold cross validation within each fold. 
CAUTION: DO NOT USE PARALLEL VERSION OF mycvFromFolds within this code
"""
from dbdKernels2 import dbdKernel, getKernelData
from myKernelData import *
from BISEPutils import *
from myPyMLassess import *
from PyML import *
from PyML.evaluators.assess import trainTest
from mpi4py import MPI
import numpy as np
import sys
from time import time
from PyML.evaluators import roc     
from dbdscrpp3 import *

if __name__=="__main__":      
    pwKfname='dbdKernel.dbk.pkl'    
    print "Loading Kernel from file:",pwKfname   
    t1=time()
    mK=getKernelData(*dbdKernel.loader(pwKfname)) #NOTE THAT OBJECT IS BEING LOADED FROM FILE
    t2=time()
    print "Pairwise kernel loaded in",t2-t1,"s"   
    ttP,dkey=getComplexFoldList(mK)    
    mK.attachLabels(Labels(mK.labels.L))    #pyml did not work with the string labels
    print "Starting CV ... "
    Clist=[0.1,1,10,100]
    Slist=[]
    for i in range(len(ttP)):
        print "##### PROCESSING FOLD #### ",i
        rmax=0
        bcc=0
        for icc,cc in enumerate(Clist):
            s=SVM()
            s.C=cc
            ttPi=ttP[0:i]
            ttPi.extend(ttP[i+1:])
            r=mycvFromFolds(s,mK,testingPatterns=ttPi)
            rauc=r.getROC()
            if rauc>rmax:
                rmax=rauc
                bcc=cc            
            #pdb.set_trace()
        s=SVM()
        s.C=bcc
        print "##### Optimal C = ",bcc
        Slist.append(s)
    r=mycvFromFolds(Slist,mK,testingPatterns=ttP)
    if r is not None:
        output = open('result.res.pkl', 'wb')
        cPickle.dump((r,dkey), output)
        output.close()
        (auc,_,_)=getAUC((r,dkey))
        print "CV Complete. Total Time taken (s):",time()-t1
        print "Complex-wise average AUC",auc
        print "Overall AUC = ", r.getROC()
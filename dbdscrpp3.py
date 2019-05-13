# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 23:32:18 2012

@author: fayyaz
Description: Performs cross validation on a kernel file
"""
from dbdKernels import dbdKernel, getKernelData
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
import random
def getComplexFoldList(mK,nfolds=4,shuffle=False):    
    """
    If nfolds is a list of lists containing the complex ids for each fold then that is used
    otherwise it is created
    """
    dkey={}    
    for i,e in enumerate(mK.labels.patternID):
        try:
            v=dkey[e[0]]
        except:
            dkey[e[0]]=([],[])            
        dkey[e[0]][0].append(i)
        dkey[e[0]][1].append(e)
    #pdb.set_trace()
    
    if type(nfolds) is type([]):
        fcids=nfolds
        nfolds=len(fcids)
        #the fcids and dkey must have exactly the same elements
        vkeys=list(set(itertools.chain(*fcids)).intersection(dkey.keys()))
        dkey2={}
        for k in vkeys:
            dkey2[k]=dkey[k]
        dkey=dkey2
    else:
        cids=dkey.keys()
        if shuffle:
            random.shuffle(cids)
        fcids=list(chunks(cids,int(np.ceil(len(cids)/float(nfolds)))))
    print fcids
    #pdb.set_trace()
    F=[[] for f in range(nfolds)]
    for f in range(nfolds):
        for cid in fcids[f]:
            if cid in dkey:
                F[f].extend(dkey[cid][0])
            else:
                print "Warning: Complex",cid,"Not found in the kernel. Will not be used in the analysis"
    return F,dkey

def getAUC(s):
    if type(s)==type(''):
        (r,dkey)=cPickle.load(open(s, "rb" ) )
    else:
        (r,dkey)=s
 
    patid=combineList(r.getPatternID())
    vkey=dict(zip(patid,range(len(patid))))
    decfn=combineList(r.getDecisionFunction())
    lblid=combineList(r.getGivenLabels())
    cids=dkey.keys()
    D=[[] for i in cids]
    L=[[] for i in cids]
    A=[[] for i in cids]
    try:
        R=getRMSDDict('shandar_rmsd.txt')
    except:
        R=None
    Rx=[[] for i in cids]
    for i,cid in enumerate(cids):
        cidx=dkey[cid]        
        if type(cidx) is tuple: #backward compatability to old results objects 
            cidx=cidx[0]
        for e in cidx:
            n=vkey[e]
            D[i].append(decfn[n])
            L[i].append(lblid[n])
        (_,_,a)=roc.roc(D[i],L[i])
        A[i]=a
        if R is not None:
            Rx[i]=R[cid]        
    (fp,tp,auc)=roc.roc_VA(zip(D,L))
    return (auc,(fp,tp),(A,Rx,D,L,cids,r,dkey))

def getRMSDDict(fname):
    D={}
    for l in open(fname,'r'):
        x=l.split()
        D[x[0][:4]]=float(x[1])
    return D
    
if __name__=="__main__":      
    pwKfname='dbdKernel_dbd3_str_tppk.dbk.pkl'
    print pwKfname
    s=SVM()
    s.C=10    
    print "Loading Kernel from file..."
    t1=time()
    mK=getKernelData(*dbdKernel.loader(pwKfname)) #NOTE THAT OBJECT IS BEING LOADED FROM FILE
    #pdb.set_trace()
    t2=time()
    print "Pairwise kernel loaded in",t2-t1,"s"   
    #pdb.set_trace()
    ttP,dkey=getComplexFoldList(mK)
    mK.attachLabels(Labels(mK.labels.L))    #pyml did not work with the string labels
    print "Starting CV ... "
    r=mycvFromFolds(s,mK,testingPatterns=ttP)
    if r is not None:
        output = open('result.res.pkl', 'wb')
        cPickle.dump((r,dkey), output)
        output.close()
        (auc,_,_)=getAUC((r,dkey))
        print "CV Complete. Total Time taken (s):",time()-t1
        print "Complex-wise average AUC",auc
        print "Overall AUC = ", r.getROC()
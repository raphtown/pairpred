# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 14:37:02 2012
Modified: Dec 31, 2012
@author: Fayyaz
Description: Code for performing cross validation analysis using our kernel in 
parallel.
"""
from PyML.evaluators.assess import trainTest
from mpi4py import MPI
import numpy as np
import sys
import itertools
import pdb
def chunks(l, n):
    """ Yield successive n-sized chunks from list l, returns list of lists."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def mycvFromFolds(classifier, data, testingPatterns, trainingPatterns=None,ofname=None ,comm=None,myid=0,nprocs=1,**args) :
    
    """perform cross validation from folds in parallel using MPI

    :Parameters:
      - `classifier` - a classifier template. 
          Can be a list of classifiers. In that case each fold is tested using its
          own classifier.
      - `data` - a dataset
      - `testingPatterns` - a list providing the testing examples for each fold
      - `trainingPatterns` - a list providing the training examples for each fold      
          optional: When not specified, the training patterns of a fold are
          the patterns not in the test set of that fold          
      - 'ofname' (optional, default = cvResults.cvr), name of the file to 
              which the results object is saved. The object can be loaded
              using loadResults(ofname)
                  from PyML.evaluators.resultsObjects import loadResults
     - comm,myid,nprocs: mpi4py parameters for performing parallel execution
         optional, Default (None,0,1). When not specified, a single process is used.
             comm = MPI.COMM_WORLD
             myid = comm.Get_rank()
             nprocs = comm.Get_size()
         
    :Keywords:
      - `intermediateFile` - a file name to save intermediate results under
        if this argument is not given, no intermediate results are saved

    :Returns:
      None for processors with ids > 1
      a results object for processor 0
    """
    ncv=type(classifier)==type([]) #pass in a list of classifiers -- Caution: DO NOT USE PARALLEL VER.
    N=len(testingPatterns) #Number of folds
    if ncv and comm is not None:
        raise Exception('Current Version of myPyMLassess does not support parallel execution when a classifier list is specified.')
    if ncv:
        assert len(classifier)==N    
    if trainingPatterns is None:
        #L=set(data.labels.patternID) #get all the labels
        L=set(itertools.chain(*testingPatterns)) #get all the labels in the training
        trainingPatterns=[]
        for f in range(N):
            trainingPatterns.append(list(L.difference(testingPatterns[f])))
        
    assert len(trainingPatterns) == len(testingPatterns) 
    args['stats'] = False
    
    if nprocs>N:
        if myid==0:
            print "Using only",N," processors out of ", nprocs
        nprocs=N
    if myid>=N:
        return None #sys.exit(0)
    csize=int(np.ceil(N/float(nprocs)))
    mytestingPatterns=list(chunks(testingPatterns,csize))[myid]
    mytrainingPatterns=list(chunks(trainingPatterns,csize))[myid]
    
    s=[[] for p in range(len(mytestingPatterns))]
    
    for idx in range(len(mytestingPatterns)):
        if ncv:
            ss=classifier[idx]
        else:
            ss=classifier
        s[idx]=trainTest(ss, data,
                                   mytrainingPatterns[idx], mytestingPatterns[idx], **args)    
    if myid!=0:
        comm.send(s, dest=0)
        return None
    if(myid==0):
        cvResults = ss.resultsObject()
        for idx in range(len(s)):
            cvResults.extend(s[idx])
        for p in range(1,nprocs):
            s=comm.recv(source=p)
            for idx in range(len(s)):
                cvResults.extend(s[idx])
        cvResults.computeStats()
        if ofname is not None:               
            print 'Result File Saved: ', ofname
            cvResults.save(ofname,'long')
        return cvResults


    

if __name__=="__main__":        
    
    from myKernelData import *
    from scipy.spatial.distance import cdist
    N=2000
    folds=5
    X=2*np.random.random((N,2))-1 #DIFFERENT PROCESSES WILL GENERATE DIFFERENT DATA
    w=np.array([1,3])
    y=2*(np.dot(w,X.T**2)>0.5)-1    
    K=np.exp(-1*cdist( X, X )**2)#np.dot(X,X.T)
    kdata = myKernelData(K)
    kdata.attachLabels(Labels(y))
    s=SVM()
    s.C=10
    r=mycvFromFolds(s,kdata,testingPatterns=list(chunks(range(N),N/folds)))    
    if r is not None:        
        s=[]
        L=r.getGivenLabels()
        F=r.getDecisionFunction()
        #r.getPatternID() contains the list
        a=zip(F,L)
        from PyML.evaluators import roc
        print roc.roc_VA(a)
        r.plotROC()


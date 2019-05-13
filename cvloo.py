# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 12:03:18 2013

@author: root
"""
from testSingleComplex_par import *
if __name__=="__main__":
    evalROC=True #Whether we want to calculate the ROC or not 
    selectAll=True  #Whether to use all examples in the complex or not
    pdbpklpath='/s/chopin/b/grad/minhas/Desktop/DBD4N/PDBPKL4'#'/media/sf_Desktop/PDBPKLP3' # 
    pwKfname='dbdKernel_seq_tppk_d3.dbk.pkl'
    odir='./DBD4LOOCV/'
    comm = MPI.COMM_WORLD
    myid = comm.Get_rank()
    nprocs = comm.Get_size()
    (_,_,E)=cPickle.load(open(pwKfname, "rb" ))
    cids=E.Pex.keys()[90:]
    if myid==0:
        print 'Complexes being processed by me are:',cids
    for (i,cid) in enumerate(cids):        
        lname=cid+'_l'
        rname=cid+'_r'
        fname=odir+cid+'.pairpred.txt'
        if myid==0:
            print 'Currently Processing: ',i,cid
        try:       
            testComplex(cid,lname,rname,evalROC,selectAll,pwKfname,pdbpklpath,comm,myid,nprocs,fname)
        except Exception as e:
            if myid==0:
                print '###PROCESSSING FAILED FOR ',i,cid,e
        comm.Barrier()
                
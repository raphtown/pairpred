# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 18:19:12 2012

@author: Fayyaz
Description: Given a kernel file (.dbK.pkl) and a complex id (cid) and the path 
to the .pdb.pkl files for the data used to construct the kernel file and the
.pdb.pkl files for the complex, this code evaluates the output for that complex id.
"""
import os
from dbdKernels import *
from PyML import *
from PyML.evaluators.assess import trainTest
from PyML.evaluators import roc
from myPyMLassess import *
from mpi4py import MPI
from time import time
from copy import deepcopy

def write2File(L,R,Scx,rZ,pwKfname,tx,auc=None,fname=None,lname='L',rname='R'):
    """
    Write the results to a file
    Format: (Example)
        AUC = 0.614399446175
        L_A R 258 68 R_A D 129 43 -0.307060015218 1
        L_A R 261 71 R_A D 129 43 -0.50449226754 1
        ...
        L_chainid <> res_name <> residue id in pdb file <> residue number in sequence or FASTA file <> residue number in L.R <>\
        R_chainid <> res_name <> residue id in pdb file <> residue number in sequence or FASTA file <> residue number in R.R <>\
        decision value <> label 
    AUC is only printed if evalROC was set to true
    labels are true labels if evalROC is true (otherwise all are -1)
    """
    lresi=dict (zip(L.resi.values(),L.resi.keys()))
    lR2Si=dict (zip(L.S2Ri,range(len(L.S2Ri))))
    rresi=dict (zip(R.resi.values(),R.resi.keys()))
    rR2Si=dict (zip(R.S2Ri,range(len(R.S2Ri))))
    lname=lname+':'
    rname=rname+':'
    if fname is None:
        fname=Scx[0][0]+'.InterPRed.txt'
    f = open(fname, 'w+')
    try:
        s='# Complex Name: '+Scx[0][0];f.write(s+'\n')        
        s='# ligand file: '+L.ifname;f.write(s+'\n')
        s='# Receptor file: '+R.ifname;f.write(s+'\n')
        s='# Ligand Name: '+lname;f.write(s+'\n')       
        s='# Receptor Name: '+rname;f.write(s+'\n')       
        s='# Time taken (s): '+str(tx);f.write(s+'\n')
        s='# Number of examples: '+str(len(Scx));f.write(s+'\n')
        if auc is not None:
            s='AUC = '+str(auc)+'\n'
            f.write(s)
        for (k,(cid,(rl,rr),rlbl)) in enumerate(Scx):
            rla=rra='-'
            rlai=rrai=-1
            if rl in lR2Si:
                rla=L.seq[lR2Si[rl]]
                rlai=lR2Si[rl]
            if rr in rR2Si:
                rra=R.seq[rR2Si[rr]]
                rrai=rR2Si[rr]
            s=lname+lresi[rl][0]+' '+rla+' '+lresi[rl][1]+' '+str(rlai)+' '+str(rl)+' '+\
             rname+rresi[rr][0]+' '+rra+' '+rresi[rr][1]+' '+str(rrai)+' '+str(rr)+\
             ' '+str(rZ[k])+' '+str(rlbl)+'\n'
            f.write(s)
        
    except Exception as e:
        print "Error saving Results: ",e
        raise e
    finally:
        f.close()


def testComplex(cid,lname,rname,evalROC,selectAll,pwKfname,pdbpklpath,comm,myid,nprocs,ofile=None):
    (Ktr,S,dK,E)=dbdKernel.loader(pwKfname)  
    if evalROC and cid not in E.Pex:
        emsg="Cannot evaluate the ROC as positive and negative examples for the specified complex are not known"
        print "Error:",emsg
        raise Exception(emsg)   
    #Construct the complex kernel
    cids=list(np.sort(E.Pex.keys())) #sending using MPI or picking changes the order of the dictionary   
    if myid==0:
        t0=time()
    cKtrtt=computecKdict(cids,pdbpklpath,E,scid=cid,comm=comm,myid=myid,nprocs=nprocs) #cid does not have to be in E
    
    if myid==0:  
        print "Time taken in evaluating Complex-Kernel (s):",time()-t0
        lpdbfname=os.path.join(pdbpklpath,cid+'_l_u.pdb.pkl')
        rpdbfname=os.path.join(pdbpklpath,cid+'_r_u.pdb.pkl')
        L=myPDB.loader(lpdbfname)
        R=myPDB.loader(rpdbfname)
        #get examples
           
        if not evalROC:
            cide=list(itertools.product(L.S2Ri,R.S2Ri)) #this way we do not get the nans
            if not selectAll:
                cide=cide[:5000]
            Sc=zip([cid]*len(cide),cide,[-1]*len(cide))            
        else:            #if we know the positive examples
            print "Original number of examples: ",len(E.Pex[cid][0])+len(E.Nex[cid])     
            cidep=E.Pex[cid][0]
            if selectAll:
                E.selectNegEx(cname=[cid],p=100.0,M=None) #select all negative examples
            ciden=E.Nex[cid]        
            Scp=zip([cid]*len(cidep),cidep,[1]*len(cidep))   
            Scn=zip([cid]*len(ciden),ciden,[-1]*len(ciden))
            Sc=Scp
            Sc.extend(Scn)            
        Scx=Sc
        S0=deepcopy(S)
        Scx_n=[]
        rZ=[]
        rL=[]
        cKtt=computecKdict([],pdbpklpath,E,scid=cid) #E is not needed at all
        print "Total number of examples: ",len(Scx)
        
        t0=time()
        #Scx,S0,cKtrtt,cKtt,dK,cid,Ktr
    else:
        Scx=None
        S0=None
        cKtrtt=None
        cKtt=None        
    (Scx,S0,cKtrtt,cKtt)=comm.bcast((Scx,S0,cKtrtt,cKtt),root=0)
    csize=int(np.ceil(len(Scx)/float(nprocs)))
    
    mScx=list(chunks(Scx,csize))[myid]
    print 'My id is:',myid,'and I am processing',len(mScx),'test examples'
    mScx_n=[]
    mrZ=[]
    mrL=[]
    if myid==0:
        t0=time()
    for idxx,Sc in enumerate(chunks(mScx,3000)):
        S=deepcopy(S0)
        #compute the tr-tt kernel
        Ktrtt,S,Sc=dbdK2pwK(cKtrtt,Sr=S,Sc=Sc)
        #print "Time taken in evaluating PW-Kernel (s):",time()-t0
        #remove nans and normalize            
        (Kttd,Sc,_)=dbdK2pwK(cKtt,Sr=Sc,justdiag=True)
        #print "NORMS=",np.linalg.norm(Ktrtt),np.linalg.norm(Kttd)
        nanidx=~np.isnan(Kttd)        
        Ktrtt=Ktrtt[:,nanidx]
        Kttd=Kttd[nanidx]
        Sc=[s for i,s in zip(nanidx,Sc)  if i]
        Ktrttn=Ktrtt/np.sqrt(dK[:,np.newaxis]*Kttd[np.newaxis,:]) #Normalization      
        Ktrtt=None
        #construct labels and kernel
        #print "NORMS = ",np.linalg.norm(Ktr),np.linalg.norm(Ktrtt),np.linalg.norm(Kttd)
        K=np.vstack((np.hstack((Ktr,Ktrttn)),np.hstack((Ktrttn.T,np.eye(Ktrttn.shape[1])))))
        Ktrttn=None;  
        trpid=[i for (i,(c,_,_)) in enumerate(S) if c!=cid] #training indices (test complex is excluded)
        S.extend(Sc)
        mK=getKernelData(K,S)
        K=None;
        #Train and test
        mK.attachLabels(Labels(mK.labels.L))
        s=SVM();s.C=10
        r=mycvFromFolds(s,mK,trainingPatterns=[trpid],testingPatterns=[range(Ktr.shape[0],len(S))])
        mK=None	
        mrZ.extend(r.getDecisionFunction())
        mrL.extend([l for (_,_,l) in Sc])
        mScx_n.extend(Sc)
            #nexproc=nexproc+len(Sc)
            #print "Classified",nexproc,"/",len(Scx),"Examples. Time:",time()-t0
    if myid==0:
        
        Scx_n=mScx_n
        rZ=mrZ
        rL=mrL
        for p in range(1,nprocs):
            (mrZ,mrL,mScx_n)=comm.recv(source=p)
            rZ.extend(mrZ)
            rL.extend(mrL)
            Scx_n.extend(mScx_n)
        tx=time()-t0
        print "Classification complete in (s):",tx
        auc=None
        if evalROC:
            (fp,tp,auc)=roc.roc(rZ,rL)
            print 'AUC=',auc
        write2File(L,R,Scx_n,rZ,pwKfname,tx,auc=auc,lname=lname,rname=rname,fname=ofile)
    else:
        comm.send((mrZ,mrL,mScx_n), dest=0)
        
if __name__=="__main__":
    lname='ISG15' #Name of the ligand protein
    rname='NS1B' #Name of the Receptor Protein
    cid='MD7B' #Complex id
    evalROC=False #Whether we want to calculate the ROC or not 
    selectAll=True  #Whether to use all examples in the complex or not
    pdbpklpath='/s/chopin/b/grad/minhas/Desktop/DBD4N/PDBPKL4'#'/media/sf_Desktop/PDBPKLP3' # 
    pwKfname='dbdKernel_dbd34_str.dbk.pkl'
    comm = MPI.COMM_WORLD
    myid = comm.Get_rank()
    nprocs = comm.Get_size()
    testComplex(cid,lname,rname,evalROC,selectAll,pwKfname,pdbpklpath,comm,myid,nprocs)

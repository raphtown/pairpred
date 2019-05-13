# -*- coding: utf-8 -*-
"""
Updated: 16JAN13 Created on Wed Oct 24 00:30:57 2012
@author: Fayyaz Minhas
Updated to execute in matrix mode (fast!)
This module implementes different kernels for the dbd data. It allows constr-
uction of pairwise kernel given the just the pdbpkl files and the getExamplesDBD
object corresponding to them.
Now it saves the kernel in two separate files (loading is backward compat.)
If you run the module then it will ecaluate the dbdKernel and save it
    mpirun -np n python dbdKernels.py
When ussing the saved dbdKernel file 
    from dbdKernels import pwKernel,dbdKernel,complexKernel
"""
from BISEPutils import *
from myPDB import *
from getExamplesDBD import *
from protKernel import *
from myKernelData import *
import cPickle
from mpi4py import MPI
from copy import deepcopy
from time import time

def computecKdict(KK,pdbpklpath,E,scid=None,comm=None,myid=0,nprocs=1):    
    """
    Computes the dictionary object cKdict (complex-wise kernel)
    Input: 
        KK: list of complex ids for which the kernel dict is to be computed
        pdbpklpath: directory containing the PDBPKL files
        E: The getExamplesDBD object, containing the examples
        scid: When not None, the complex kernel is evaluated between complexes in KK and scid only
        (comm,myid,nprocs) mpi4py parallelization stuff
    Returns:
        cKdict: A dictionary object with (cidA,cidB) as the key and the value is the complexKernel object between these complexes        
            
    """
    aA=False 
    if len(KK)==0:
        assert scid is not None
        KK=[scid]
        aA=True
    Nc=len(KK)     
    bA=False
    if scid is None:
        clist=[(a,b) for a in range(Nc) for b in range(a,Nc)]
    else:
        clist=[(a,scid) for a in range(Nc)]
        bA=True
    
    csize=int(np.ceil(len(clist)/float(nprocs)))
    gclist=list(chunks(clist,csize))
    myclist=gclist[myid]
    mycK={}
    for (cpa,cpb) in myclist:
        if bA:
            kcpb=cpb
        else:
            kcpb=KK[cpb]
        cpkey=(KK[cpa],kcpb)
        print "Processing complexes : ",cpkey
        mycK[cpkey]=complexKernel(KK[cpa],kcpb,pdbpklpath,E,aAll=aA,bAll=bA)
    
    cKdict=None  
    if(myid!=0):
        comm.send(mycK, dest=0)
        
    if(myid==0):
        gcK=[mycK]
        for p in range(1,nprocs):
            gcK.append(comm.recv(source=p))            
        #self.E=E
        #self.pdbpklpath=pdbpklpath
        cKdict=mergeDicts(gcK)
        
    return cKdict

def S2D(Sr,cKdict,first=True):
    """
    Convers the example list to a dictionary.
    Inputs:
        Sr: Input example list with each entry being a tuple (cid,(lrid,rrid),lbl)
            where   cid: Complex id
                    lrid: ligand (L) residue id in "L.R"
                    rrid: receptor (R) residue id  in "R.R"
                    lbl: label +1/-1
        cKdict: complex kernel dictionary object
        first: (optional,default "True")
            True: The indices from the first complex in the complex-wise kernel are used
            False: The indices from the second complex in the complex-wise kernel are used (required when Sc is specified in dbdK2pwK)
                
    Returns:
        dc: dictionary object dc[cid] which contains cid four lists as follows:
            dc[cid][0]: indices of residues for the ligand in the complex kernel for cid (note: indices are in the complex kernel and not in L.R)
            dc[cid][1]: indices of residues for the receptor in the complex kernel for cid (note: indices are in the complex kernel and not in R.R)
            dc[cid][2]: labels of examples 
            dc[cid][3]: indices of examples in the pairwise kernel
        The example (cid,(dc[cid][0][i],dc[cid][1][i]),dc[cid][2][i]) is mapped to the pairwise kernel at index dc[cid][3][i]
        Srn: same as Sr except that the examples are lined up such that example Srn[i] corresponds to the index i in pairwise kernel K[i,:]
    """
    cids=list(set([cid for (cid,(r,c),lbl) in Sr]))
    cids.sort() ##VERY IMPORTANT
    dc={}
    l0=0
    Srn=[]   
    for cid in cids:        
        if cid not in dc:
            dc[cid]=[[],[],[],[]]
        if first:           
            if (cid,cid) in cKdict:
                tK=cKdict[(cid,cid)]            
            else:
                c2=[c2 for (c1,c2) in cKdict.keys() if c1==cid]                
                tK=cKdict[(cid,c2[0])]
            dc[cid][0]=[tK.lAidx[l] for (x,(l,r),lbl) in Sr if x==cid]
            dc[cid][1]=[tK.rAidx[r] for (x,(l,r),lbl) in Sr if x==cid]
        else:
            c1=[c1 for (c1,c2) in cKdict.keys() if c2==cid]
            tK=cKdict[(c1[0],cid)] 
            dc[cid][0]=[tK.lBidx[l] for (x,(l,r),lbl) in Sr if x==cid]
            dc[cid][1]=[tK.rBidx[r] for (x,(l,r),lbl) in Sr if x==cid]            
        dc[cid][2]=[lbl for (x,(l,r),lbl) in Sr if x==cid]
        Srn.extend([(x,(l,r),lbl) for (x,(l,r),lbl) in Sr if x==cid])
        l1=l0+len(dc[cid][2])
        dc[cid][3]=range(l0,l1)
        l0=l1                
    return (dc,Srn)
    
def dbdK2pwK(cKdict,Sr,Sc=None,justdiag=False):
    """
    Converts a dictionary of complex kernels to a pairwise kernel
    Input:
            cKdict: obtained using computecKdict
            Sr: list of Examples for the rows with each entry being a tuple (cid,(lrid,rrid),lbl)
            Sc: list of examples for the cols (when None (default): taken to be same as Sr)
            justdiag: set to 'True' If only the diagonal is to be computed
            
    Return:
            Krc: The pairwise kernel
            Sr: list of examples for the rows 
            Sc: list of examples for the cols
            The complexes are alphabetically sorted in both Sr and Sc so the order of examples is changed from the input
        Krc[i,j] represents the kernel value between examples Sr[i] and Sc[j]
        When justDiag is True a 1D array of the diagonal values is returned
    """
    def getPWK(_cida,_cidb):
        """
        evaluates the pairwise kernel (TPPK, MLPK and SUM) between two complexes _cida,_cidb
        Requires that the variable tK contains the complex kernels between the complexes and
        that the examples appear in dictionaries dr and dc. _cida is taken to be the rows
        and _cidb is taken to be the columns.
        """
        krl=tK.K_rl[dr[_cida][1],:][:,dc[_cidb][0]]
        krr=tK.K_rr[dr[_cida][1],:][:,dc[_cidb][1]]
        kll=tK.K_ll[dr[_cida][0],:][:,dc[_cidb][0]]
        klr=tK.K_lr[dr[_cida][0],:][:,dc[_cidb][1]]
        #k_o=krl+krr+kll+klr
        k_tppk=kll*krr+klr*krl
        #k_mlpk=(kll+krr-klr-krl)**2
        kv=k_tppk
        return kv    
    (dr,Sr)=S2D(Sr,cKdict)   
    sym=False
    if Sc is None:
        Sc=Sr
        dc=dr
        sym=True
    else:
        (dc,Sc)=S2D(Sc,cKdict,first=False) 
    if justdiag:
        assert sym
        Krc=np.zeros(len(Sr),np.single)    
    else:
        Krc=np.zeros((len(Sr),len(Sc)),np.single)
    st=0    
    for aidx,cida in enumerate(dr.keys()):    
        if justdiag:
            tK=cKdict[(cida,cida)]
            kab=getPWK(cida,cida).diagonal()
            #pdb.set_trace()
            Krc[dr[cida][3]]=np.single(kab)
            continue
        if sym:
            st=aidx
        for cidb in dc.keys()[st:]:    
            #print "PW Processing",cida,cidb                   
            try:
                tK=cKdict[(cida,cidb)]
                _cida,_cidb=cida,cidb
            except Exception as exc:
                if sym:
                    _cidb,_cida=cida,cidb
                    tK=cKdict[(_cida,_cidb)]
                else:
                    raise exc            
            kab=getPWK(_cida,_cidb)
            if _cida!=cida:
                kab=kab.T
            Krc[np.ix_(dr[cida][3],dc[cidb][3])]=kab
            if sym:
                Krc[np.ix_(dc[cidb][3],dr[cida][3])]=kab.T                  
    return (Krc,Sr,Sc)       

class dbdKernel:
    """
    Given a getExamplesDBD object, it evaluates the kernel values over all
    examples from all complexes in it
    E: the getExamplesDBD object
    pdbpklpath: path where the pdbpkl files are kept
    ofname: where the output file is stored (default: dbdKernel.dbK.pkl)
    cKdict: is a dictionary object which when given the (cidA,cidB) key will 
        return the complex kernel between them
        Since the kernel is evaluated between every two complexes so only the
        upper triangular matrix (along with the diagnoal) is stored
        so if you cannot find (cidA,cidB) try (cidB,cidA)
    cids: complex ids for which the kernel has been evaluated (obtained from E)
    Note: This object support MPI style parallelization and all you have to 
        do is to call the constructor in a program that has been run using 
        mpiexec or mpirun. Only processor 0 will contain the 'valid' (complete)
        object which is saved in ofname. Each process gets equal load.
    """
    def __init__(self,pdbpklpath,E0,ofname=None,incids=None):
        comm = MPI.COMM_WORLD
        myid = comm.Get_rank()
        nprocs = comm.Get_size()        
        E=getValidE(E0,pdbpklpath,incids) 
        self.E=E
        KK=E.Pex.keys()
        if myid==0:
            print "Number of Examples (+,-)",np.sum([len(E.Pex[k][0]) for k in E.Pex]),np.sum([len(E.Nex[k]) for k in E.Nex])
            txx0=time()
            ignored=list(set(E0.Pex.keys()).difference(KK))
            if(len(ignored)):
                print 'Ignoring the following complexes because at least one of the \
                four required pdbpkl files could not be found in', pdbpklpath,':'\
                ,ignored
            print "Evaluating kernel over: ", KK
            print "Using",nprocs,"processe(s)."
        #pdb.set_trace()
        cKdict=computecKdict(KK,pdbpklpath,E,comm=comm,myid=myid,nprocs=nprocs)
        if myid==0:            
            self.__dbdK2pwK__(KK,E,cKdict)
            print "TIME TAKEN FOR KERNEL EVALUATION (s): ",time()-txx0
            txx0=time()
            if ofname is None:
                ofname='dbdKernel.dbK.pkl'            
            #self.ofname=ofname         
            self.save(ofname)
            print "TIME TAKEN FOR Saving kernel file (s): ",time()-txx0
            

            
    def __dbdK2pwK__(self,cids,E,cKdict):        
        S=[(cid,e,+1) for cid in cids for e in E.Pex[cid][0]]
        N=[(cid,e,-1) for cid in cids for e in E.Nex[cid]]
        S.extend(N)
        t0=time()
        (K,S,_)=dbdK2pwK(cKdict,S)    
        print "PW KERNEL TIME: ",time()-t0
        #HANDLE NANS HERE: If a diagonal is a nan, remove it
        nanidx=~np.isnan(K.diagonal())
        K=K[nanidx,:]
        K=K[:,nanidx]
        S=[s for (s,n) in zip(S,nanidx) if n]
        #Normalize
        dK=np.diag(K)
        K=K/np.sqrt(dK[:,np.newaxis]*dK[:,np.newaxis].T)
        self.K=K
        self.S=S
        self.dK=dK


    def save(self,ofname=None):
        """
        Save the object
        """
        if ofname is None:
            ofname='dbdKernel.dbK.pkl'
        output = open(ofname, 'wb')
        cPickle.dump((self.S,self.dK,self.E), output,-1)        
        output.close()
        output = open(ofname+'.np', 'wb')
        np.save(output,self.K)
        output.close()
        print "File Saved: "+ofname
    @classmethod  
    def loader(self,pklfile):
        """
        Load the class object from a pickel file
        """        
        z=cPickle.load(open(pklfile, "rb" ))
        if len(z)==3: #FOR BACKWAED COMPATABILITY
            K=np.load(pklfile+'.np')
            z=(K,)+z
        return z
        
def getKernelData(K,S,dK=None,E=None):
    """
    Return the myKernelData object (convert to pyml style kernel)
    """
    #pdb.set_trace()

    mK=myKernelData(np.double(K),patternID=zip(*zip(*S)[:2]),labels=list(zip(*S)[-1]))
    return mK
        
class complexKernel:
    """
    Evaluate the kernel values between 2 complexes A & B as follows
    (let l and r represent the ligand and receptor within each complex)
        K_ll=K(lA,lB)
        K_lr=K(lA,rB)
        K_rl=K(rA,lB)
        K_rr=K(rA,rB)
    It evaluates these kernels only over examples passed in E
    lAidx,rAidx,lBidx,rBidx are dictionary objects which when given the 
    residue index will return the index of that residue in the respective kernel
    Note: the inner kernel is computed using protKernel
    """
    def __init__(self,cidA,cidB,pdbpklpath,E,aAll=0,bAll=0):
        self.cidA=cidA
        self.cidB=cidB
        lAname=os.path.join(pdbpklpath,cidA+'_l_u.pdb.pkl')
        rAname=os.path.join(pdbpklpath,cidA+'_r_u.pdb.pkl')
        lBname=os.path.join(pdbpklpath,cidB+'_l_u.pdb.pkl')
        rBname=os.path.join(pdbpklpath,cidB+'_r_u.pdb.pkl')        
        lA=myPDB.loader(lAname)
        rA=myPDB.loader(rAname)
        lB=myPDB.loader(lBname)
        rB=myPDB.loader(rBname)
        #pdb.set_trace()
        if aAll:
            lAidx=range(len(lA.R))
            rAidx=range(len(rA.R))
        else:
            lAidx_p,rAidx_p=tuple(map(np.unique,zip(*E.Pex[cidA][0])))
            lAidx_n,rAidx_n=tuple(map(np.unique,zip(*E.Nex[cidA])))
            lAidx=list(np.unique(np.concatenate((lAidx_p,lAidx_n))))
            rAidx=list(np.unique(np.concatenate((rAidx_p,rAidx_n))))
        if bAll:            
            lBidx=range(len(lB.R))
            rBidx=range(len(rB.R))
        else:
            lBidx_p,rBidx_p=tuple(map(np.unique,zip(*E.Pex[cidB][0])))
            lBidx_n,rBidx_n=tuple(map(np.unique,zip(*E.Nex[cidB])))
            lBidx=list(np.unique(np.concatenate((lBidx_p,lBidx_n))))
            rBidx=list(np.unique(np.concatenate((rBidx_p,rBidx_n))))
        self.K_ll=protKernel(lA,lB,lAidx,lBidx).K
        self.K_lr=protKernel(lA,rB,lAidx,rBidx).K
        self.K_rl=protKernel(rA,lB,rAidx,lBidx).K
        self.K_rr=protKernel(rA,rB,rAidx,rBidx).K
        print "#NANS in kernels: ",np.sum(np.isnan(self.K_ll)),np.sum(np.isnan(self.K_lr)),np.sum(np.isnan(self.K_rl)),np.sum(np.isnan(self.K_rr))
        self.lAidx=dict(zip(lAidx,range(len(lAidx))))#lAidx;#
        self.rAidx=dict(zip(rAidx,range(len(rAidx))))#rAidx;#
        self.lBidx=dict(zip(lBidx,range(len(lBidx))))#lBidx;#
        self.rBidx=dict(zip(rBidx,range(len(rBidx))))#rBidx;#
    def __str__(self):
        return "complexKernel Object over complexes :"+self.cidA+','+self.cidB
        
def getValidE(E0,pdbpklpath,cids=None):
    """
    All 4 files for the two complexes must exist for it to be included in the 
    training. This function removes all those complexes from E0 which do not satisfy this 
    requirement.
    """
    E=deepcopy(E0)
    K=set(E.Pex.keys())
    import os
    import glob
    s=set()
    for k in K:
        f=os.path.join(pdbpklpath,k+"*.pdb.pkl")
        if len(glob.glob(f))==4:
            s.add(k)
    if cids is not None:
        s=s.intersection(set(cids))
    tlist=list(s)
    E.Pex=copy_dict(E.Pex,*tlist)
    E.Nex=copy_dict(E.Nex,*tlist)
    return E
    
if __name__=="__main__":    
    pdbpklpath='/s/chopin/b/grad/minhas/Desktop/DBD4N/PDBPKL4'#'./DBD3PKL'#'/s/chopin/b/grad/minhas/PDBPKLP3'#'/s/chopin/b/grad/minhas/Desktop/PDBPKLA' #
    exfname=pdbpklpath+'/E_6.0.lbl.pkl' #    '/media/sf_Desktop/PIANO/PDBPKL/E_6.0.lbl.pkl'# 
    ofname='dbdKernel.dbk.pkl'
    
    E=getExamplesDBD.loader(exfname)
    f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN', '1N2C', '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
    f4=['2A5T', '3CPH', '1ZHH', '2ABZ', '1LFD', '2OUL', '1JIW', '2B4J', '1SYX', '1FLE', '1JTG', '2AYO', '4CPA', '1CLV', '1OC0', '1XU1', '1R6Q', '2O3B', '1US7', '3D5S', '1JZD', '1HCF', '1OYV', '2OZA', '1H9D', '2A9K', '2J0T', '2Z0E', '3BP8', '2IDO', '1WDW', '1ZLI', '2VDB', '1RV6', '1FFW', '1F6M', 'BOYV', '1JWH', '2OOR', '1MQ8', '1GL1', '1PVH', '2I9B', '1OFU', '1GXD', '3SGQ', '1JK9', '1ZM4', '1FCC', '2G77', '2J7P', '2FJU']
    incids=f3+f4
    #pdb.set_trace()
#    skeys=E.Pex.keys()[0:10]
#    E.Pex=copy_dict(E.Pex,*skeys)
#    E.Nex=copy_dict(E.Nex,*skeys)
#    tlist=['1A2K','1ACB','1AHW','1XQS']
#    E.Pex=copy_dict(E.Pex,*tlist)
#    E.Nex=copy_dict(E.Nex,*tlist)
    
#    print E.Pex.keys()
    dK=dbdKernel(pdbpklpath,E,ofname,incids)
    
    #comm = MPI.COMM_WORLD
    #myid = comm.Get_rank()
    #if myid==0:  
   #     dK.save()
        #pK=pwKernel.loader(ofname) #NOTE THAT OBJECT IS BEING LOADED FROM FILE
   #     pK=pwKernel(dK)
     #   mK=pK.getKernelData()
   #     plt.imshow(mK.getKernelMatrix());plt.show()
    #cK=complexKernel('1JPS','2HMI',pdbpklpath,E)
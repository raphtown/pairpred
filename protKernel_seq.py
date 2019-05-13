# -*- coding: utf-8 -*-
"""
Created on Thu Oct 04 01:16:00 2012
@author: Afsar
"""
from BISEPutils import *
from myPDB import *
from time import time
from scipy.spatial.distance import cdist
from Bio.SeqUtils.ProtParam import ProteinAnalysis
def findProteinRASA(L):
    """
    Given a myPDB object, it computes the relative accessible surface area for
    the whole protein
    """
    a=L.ASA[L.S2Ri]
    m=ProteinAnalysis(L.seq).molecular_weight()
    easa=(3.6162*(m**0.792))
    return np.sum(a[a==a])/easa
    
    #return np.sum(L.ASA)
def findNASA(L):
    """
    Finds the accessible surface area in the neighborhood of the residue and return
    the normalized mean
    """
    Z=np.zeros(len(L.R))
    for i in range(len(L.R)):
        a=L.ASA[L.S[0][i]]
        Z[i]=(np.mean(a[a==a])-41.6)/17.6 #0 mean and 1.0 std
    return Z
def processPSAIA(L):
    """
    Normalizes all PSAIA features to 0-1
    """
    RDPX_M=np.array([7.5131,2.4013,7.658,1.8651,8.0278,7.0175])
    RDPX_m=-1.0
    RCX_M=np.array([11.153,4.8229,12.212,4.4148,19.84,9.4199])
    RCX_m=-1.0
    RRASA_M=np.array([167.0,368.04,124.66,316.96,253.85])
    RRASA_m=0.0
    CASA_m=np.array([1412.0,235.32,1174.2,362.0,940.64])
    CASA_M=np.array([36661.0,7584.9,29756.0,15550.0,21489.0])
    RASA_m=0.0
    RASA_M=np.array([273.22,134.51,216.4,173.26,185.47])
    RHPH_m=-4.5
    RHPH_M=+4.5
    #pdb.set_trace()
    rdpx=(L.psaiaf['rdpx']-RDPX_m)/(RDPX_M-RDPX_m)[:,np.newaxis]
    rcx=(L.psaiaf['rcx']-RCX_m)/(RCX_M-RCX_m)[:,np.newaxis]
    rrasa=(L.psaiaf['rrasa']-RRASA_m)/(RRASA_M-RRASA_m)[:,np.newaxis]
    casa=(L.psaiaf['casa']-CASA_m[:,np.newaxis])/(CASA_M-CASA_m)[:,np.newaxis]
    rasa=(L.psaiaf['rasa']-RASA_m)/(RASA_M-RASA_m)[:,np.newaxis]
    rhph=(L.psaiaf['rhph']-RHPH_m)/(RHPH_M-RHPH_m)
    fv=np.vstack((rcx,rhph))
    #pdb.set_trace()
    return fv#fv/(1e-10+np.sum(np.abs(fv)**2,axis=0)**(1./2))
    
    
def processPSSM(p):
    """
    return normalized pssm fv
    """
    z=p
    return z/(1e-10+np.sum(z**2,axis=0)**(0.5))
def getCombinedFV(L):
    """
        Features
            sequence L.FV.toarray() (divide by sqrt(2) to make the norm 1)
            L.RD[0] 
            L.RD[1]
            L.rASA
            L.ASA           
            L.rasa predicted 
            L.UC
            L.DC
            findNASA(L)
            findProteinRASA(L)
            pssm
            psfm
        """
    pass
    
class protKernel:
    def __init__(self,L,R,Lidx=None,Ridx=None):
        
        """
        Lidx and Ridx are indices of the residues for which the kernel values
        are to be evaluated. The inner kernels are evaluated over all the indices
        as it takes very small amount of time
        """
        self.ss=L.ifname==R.ifname
        self.L=L
        self.R=R
        self.Lidx=Lidx
        self.Ridx=Ridx
        if Lidx is None:
            Lidx=range(len(L.R))
        if Ridx is None:
            Ridx=range(len(R.R))
        #Ks0=np.dot(np.array(L.FV.todense()).T,np.array(R.FV.todense()))
        #Kprasa=np.exp(-10*(findProteinRASA(L)-findProteinRASA(R))**2)
        #Ks=np.exp(-0.5*cdist(L.FV.toarray().T,R.FV.toarray().T)**2)
        #Kss=np.dot(L.SS.T,R.SS)     
        #Kuc=np.exp(-3*cdist(L.UC.T,R.UC.T)**2)#np.dot(L.UC.T,R.UC)
        #Kdc=np.exp(-3*cdist(L.DC.T,R.DC.T)**2)#np.dot(L.DC.T,R.DC)
        #Krd=(np.exp(-0.1*cdist(L.RD[0][:,np.newaxis],R.RD[0][:,np.newaxis])**2)+\
        #    np.exp(-0.1*cdist(L.RD[1][:,np.newaxis],R.RD[1][:,np.newaxis])**2))/2.0
        #KrASA=np.exp(-3.0*cdist(L.rASA[:,np.newaxis],R.rASA[:,np.newaxis])**2)
        #Knasa=np.exp(-1.0*cdist(findNASA(L)[:,np.newaxis],findNASA(R)[:,np.newaxis])**2)
        Kpssm=np.exp(-0.5*cdist(processPSSM(L.wpssm).T,processPSSM(R.wpssm).T))
        Kpsfm=np.exp(-0.5*cdist(processPSSM(L.wpsfm).T,processPSSM(R.wpsfm).T))
        #Kinfo=L.info[:,np.newaxis]*R.info[:,np.newaxis].T
        #Kinfo=np.exp(-3*cdist(L.info[:,np.newaxis],R.info[:,np.newaxis])**2)
        #Khseaac=np.exp(-0.5*cdist(np.vstack((L.UC,L.DC)).T,np.vstack((R.UC,R.DC)).T)**2)
        #Ksurfacc=np.exp(-3*cdist(np.vstack((L.RD[0]/8.0,L.RD[1]/8.0,L.rASA)).T,np.vstack((R.RD[0]/8.0,R.RD[1]/8.0,R.rASA)).T)**2)
        #Kpsaia=np.exp(-1.0*(cdist(processPSAIA(L).T,processPSAIA(R).T)**2))
        #Kin=Kpssm+Kpsfm+Khseaac+Ksurfacc#Ks+Kpssm+Kinfo+Kuc+Kdc+Krd+KrASA+Kprasa
        Krasa=np.exp(-3.0*cdist(L.rasa[:,np.newaxis],R.rasa[:,np.newaxis])**2)
        Kin=Kpssm+Kpsfm+Krasa
        #pdb.set_trace()
        #Kin=np.multiply(Kin,Ka)
        self.K=Kin[np.ix_(Lidx,Ridx)]
        #self.K=K
        #self.Kin=Kin
        
    def save(self,ofname):
        output = open(ofname, 'wb')
        cPickle.dump(self, output)
        output.close()
    @classmethod   
    def loader(self,pklfile):
        return cPickle.load(open(pklfile, "rb" ) )
        
    def saveForViz(self,ofname):
        d=dict()
        d['K']=self.K
        d['Ks']=None
        d['Ka']=None
        d['Kss']=None
        d['Kin']=None
        d['Lresi']=self.L.resi
        d['Rresi']=self.R.resi
        output = open(ofname, 'wb')
        cPickle.dump(d, output)
        output.close()
        
    def visualize(self,show=1):   
        lname='Ltemp.pdb'
        rname='Rtemp.pdb'
        kname='Ktemp.pkl'
        self.L.toFile(lname)        
        self.R.toFile(rname)        
        self.saveForViz(kname)
        r=""        
        if not self.ss:
            r=" -r "+rname
        execStr="python visK.py -l "+lname+r+" -k "+kname        
        if show:            
            try:
                os.system(execStr)
            except:
                print "Error Launching PyMOL."
                return
        
if __name__=="__main__":
    print "hello"
    pdbpklpath='./'#'/s/chopin/b/grad/minhas/PDBPKLP3/'#'/s/chopin/b/grad/minhas/Desktop/PDBPKLA/'#''#
    #pdbid="1SBB"
    fname_L=pdbpklpath+'2B42'+"_"+"l_u.pdb.pkl"
    fname_R=pdbpklpath+'2B4s'+"_"+"l_u.pdb.pkl"    
    L=myPDB.loader(fname_L)
    R=myPDB.loader(fname_R)
    t0 = time()
    K0=protKernel(L,R) 
    t1 = time()
    print 'Kernel eval takes %f' %(t1-t0),np.sum(np.isnan(K0.K))
#Lidx=[69]
#Ridx=[83]
#K1=protKernel(L,R,Lidx,Ridx) 
#print K1.K,K0.K[Lidx,Ridx]
#
#t0 = time()
#t1 = time()
#    print 'Kernel eval takes %f' %(t1-t0)
#    print 'Max Error is:',np.max(np.abs(K1.K-K0.K[Lidx,Ridx]))
    #K0.visualize()
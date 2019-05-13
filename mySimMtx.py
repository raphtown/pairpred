# -*- coding: utf-8 -*-
"""
Created on Thu Oct 04 01:16:00 2012
@author: Afsar
"""
from BISEPutils import *
from myPDB import *
#import subprocess
from time import time

class mySimMtx:
    def __init__(self,L,R=None):        
        #D=getPWDist(L.Coords,R.Coords) # we can get the labeling from here
        self.ss=0
        if R is None:
            R=L
            self.ss=1 #evaluating the kernel for itself
           
        #Ka=L.rASA[:,np.newaxis]*R.rASA[np.newaxis,:]
        Ka=np.exp(-5*((L.rASA[:,np.newaxis]-R.rASA[np.newaxis,:])**2))
        #Ks=np.dot(L.DC.T,R.DC)#np.dot(L.DC.T,R.DC)+np.dot(L.UC.T,R.UC)#np.multiply(np.dot(L.UC.T,R.UC),np.dot(L.DC.T,R.DC))#
        #Ksd=Ks.diagonal()
        #Ks=np.divide(Ks,np.sqrt(Ksd[np.newaxis,:]*Ksd[:,np.newaxis]))
        FVL=getSeqFV(L.stx,L.R)
        FVR=getSeqFV(R.stx,R.R)
        Ks=np.array((FVL.T*FVR).todense())
        Kss=np.dot(L.SS.T,R.SS)        
        K=np.empty((len(L.S[0]),len(R.S[0])))
        K.fill(np.nan)
        Kin=np.multiply(Kss,Ka)        
        Kin=np.multiply(Kin,Ks)
        for lidx in range(len(L.S[0])):  
            kin_l=Kin[L.S[0][lidx],:]
            if self.ss:
                range_R=range(lidx,len(R.S[0]))
            else:
                range_R=range(len(R.S[0]))
            for ridx in range_R:                
                z=np.outer(L.S[1][lidx],R.S[1][ridx])               
                kin=kin_l[:,R.S[0][ridx]] 
                K[lidx,ridx]=np.multiply(kin,z).sum()/z.sum()
                if self.ss:
                    K[ridx,lidx]=K[lidx,ridx]
        #K=1+(FVL.T*FVR).todense()  
        #Kd=K.diagonal()
        #K=np.divide(K,np.sqrt(Kd[np.newaxis,:]*Kd[:,np.newaxis]))
        self.L=L
        self.R=R
        self.K=K
        self.Ks=Ks
        self.Ka=Ka
        self.Kss=Kss
        self.Kin=Kin
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
        d['Ks']=self.Ks
        d['Ka']=self.Ka
        d['Kss']=self.Kss
        d['Kin']=self.Kin
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
            #subprocess.Popen(execStr)
            try:
                os.system(execStr)
            except:
                print "Error Launching PyMOL."
                return
        
if __name__=="__main__":
    print "hello"
    pdbpath='.//benchmark4//structures//Benchmark_4_updated//'    
    pdbid="1GL1"
    fname_L=pdbpath+pdbid+"_"+"r_b.pdb"
    fname_R=pdbpath+pdbid+"_"+"r_b.pdb"
    
    L=myPDB(fname_L)
    
    if(fname_L!=fname_R):
        R=myPDB(fname_R)
        t0 = time()
        K0=mySimMtx(L,R) 
        t1 = time()
    else:    
        t0 = time()
        K0=mySimMtx(L)    
        t1 = time()
    print 'Kernel eval takes %f' %(t1-t0)
    K0.visualize()
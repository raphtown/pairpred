# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 13:16:51 2012

@author: root

Description:
    make a plot of accuracy vs. ASA
"""
from dbdscrpp3 import getAUC
from PyML.evaluators import roc
from dbdKernels2 import *
from myPDB import *
from scipy.stats.stats import pearsonr    

def getSamplePoints(X,p):
    h,b=np.histogram(X,bins=100*p)
    hc=np.zeros(len(b))
    hc[0]=0
    hc[1:]=np.cumsum(h)
    x=np.linspace(0,len(X),p)
    x0=np.interp(x, hc, b)
    return x0
def rasaPlot(rAsa,Dxx,Lxx,Np=10):
    #nb=getSamplePoints(rAsa,Np)
    nb=np.linspace(np.min(rAsa),np.max(rAsa),Np)
    R=np.zeros(len(nb)-1)
    R.fill(np.nan)
    xx=np.zeros(len(nb)-1)
    pp=np.zeros(len(nb)-1)
    nn=np.zeros(len(nb)-1)
    for idx in range(len(nb)-1):
        vidx=np.logical_and(rAsa>=nb[idx], rAsa<nb[idx+1])
        v=Dxx[vidx]
        l=Lxx[vidx]
        xx[idx]=(nb[idx]+nb[idx+1])/2.0
        pp[idx]=np.sum(l==+1)
        nn[idx]=np.sum(l!=+1)
        if(pp[idx]>10 and nn[idx]>10):
            try:
                (_,_,R[idx])=roc.roc(list(v),list(l))
            except:
                continue   
    naidx=~np.isnan(R)
    print "Correlation Coefficient: ", pearsonr(xx[naidx],R[naidx])
    #PLOTTING ONLY CODE
    stat="$\Delta$rASA"
    ww=np.diff(nb)
    plt.figure(0)
    plt.plot(xx,R,'-o')
    plt.xlabel(stat)
    plt.ylabel("AUC")
    plt.title('AUC vs $\Delta$rASA')
    plt.grid()
    plt.figure(1)
    plt.plot(nb[:-1],pp,'r-^',label="Interacting residues")
    plt.plot(nb[:-1],nn,'k-v',label="Not-Interacting residues")
#    plt.bar(nb[:-1],pp,color='r',width=ww,label="+1")
#    plt.bar(nb[:-1],nn,color='k',width=ww,bottom=pp,label="-1")
    plt.xlabel(stat)
    plt.ylabel("Counts")
    plt.legend(loc=0)
    plt.title('Number of residues vs. $\Delta$rASA')
    plt.grid()
    plt.show()
def getDASA(uname,bname):
    """
    Given the pdb.pkl file names of the bound and unbound forms of a protein,
    this code returns ASA (or any other stat) for both the proteins after 
    "aligning" the two vectors
    """
    if type(uname)==type(''):
        L=myPDB.loader(uname) 
    else:
        L=uname
    if type(bname)==type(''):
        R=myPDB.loader(bname)
    else:
        R=bname
    (u2b,b2u)=mapU2B(L.seq,L.S2Ri,len(L.R),R.seq,R.S2Ri,len(R.R))
    
    lasa=L.rASA
    xasa=R.rASA
    #lM=ProteinAnalysis(L.seq).molecular_weight()
    #rM=ProteinAnalysis(R.seq).molecular_weight()
    #pdb.set_trace()    
    
    rasa=np.zeros(len(lasa))
    rasa.fill(np.nan)
    for i in range(len(L.R)):
        if u2b[i] is not np.nan:
            rasa[i]=xasa[u2b[i]]
    return (lasa,rasa)
if __name__=="__main__":
    fname='../Results/result_86.9n.res.pkl'
    pdbpklpath='../PDBPKLP3'
    (auc,(fp,tp),(A,Rx,Dx,Lx,cids,r,dkey))=getAUC(fname)
    F=[[] for c in cids]
    for i,cid in enumerate(cids):
        print 'Processing',cid
        L=myPDB.loader(os.path.join(pdbpklpath,cid+'_l_u.pdb.pkl'))
        R=myPDB.loader(os.path.join(pdbpklpath,cid+'_r_u.pdb.pkl'))
        Lb=myPDB.loader(os.path.join(pdbpklpath,cid+'_l_b.pdb.pkl'))
        Rb=myPDB.loader(os.path.join(pdbpklpath,cid+'_r_b.pdb.pkl'))
        lurasa,lbrasa=getDASA(L,Lb)
        rurasa,rbrasa=getDASA(R,Rb)
        ldasa=np.abs(lurasa-lbrasa)
        rdasa=np.abs(rurasa-rbrasa)
        for (_,(lidx,ridx)) in dkey[cid][1]:
            f=ldasa[lidx]+rdasa[ridx]#L.rASA[lidx]+R.rASA[ridx]#len(L.S[0][lidx])+len(R.S[0][ridx])#L.psaiaf['rhph'][lidx]+R.psaiaf['rhph'][ridx]#-(L.RD[0,lidx]+L.RD[1,lidx]+R.RD[0,ridx]+R.RD[1,ridx])
            F[i].append(f)
    
    (fp,tp,auc)=roc.roc_VA(zip(F,Lx))    
    plt.plot(fp,tp);plt.xlabel('FPR');plt.ylabel('TPR');plt.title('ROC for $\Delta$rASA. AUC = '+str(auc));plt.ylim([0,1]);plt.show()
    Fx=list(itertools.chain(*F))
    Dxx=list(itertools.chain(*Dx))
    Lxx=list(itertools.chain(*Lx))
    nidx= ~np.isnan(Fx)    
    rasaPlot(np.array(Fx)[nidx],np.array(Dxx)[nidx],np.array(Lxx)[nidx],Np=10)
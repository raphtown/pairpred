# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:07:24 2013
Code for comparing the performance of LOOCV over all examples with the one 
over selected examples
@author: root
"""

from analyzePredFile import *
from BISEPUtils import getFileParts
from dbdscrpp3 import getAUC
import glob
from PyML.evaluators import roc   
from scipy.spatial.distance import cdist
import pdb  
def getAUC4Protein(lrV):
    vl=map(list, zip(*lrV.values()));vv=vl[0];ll=vl[1]    
    (_,_,a)=roc.roc(vv,ll)
    vv=np.array(vv)
    ll=np.array(ll)
    return (a,vv,ll)
def findNTPinTop(Mvx,Mlx,Mvshape,top):
    #find the number of true positives in the top 'top'
    sidx=np.argsort(Mvx)[::-1]
    L=[si for si,i in enumerate(sidx[:top]) if Mlx[i]==1]
    dntp=[] #distance from the nearest true positive
    #find the rank of the first positive
    
    (rv,cv)=np.unravel_index(sidx[:top], Mvshape)
    (rl,cl)=np.unravel_index(np.nonzero(Mlx==1), Mvshape)
    rcl=np.array((rl.flatten(),cl.flatten()))
    rcv=np.array((rv.flatten(),cv.flatten()))
    D=cdist(rcl.T, rcv.T, 'cityblock')
    di=np.argmin(D,axis=0)
    dntp=np.array([D[dix,i] for i,dix in enumerate(di)])    
    if len(L):
        fpi=L[0]        
    else:
        for si,i in enumerate(sidx[top:]):
            if Mlx[i]==1:
                break
        fpi=top+si
            
    ttp=len(L) #top true positives

    
    return (ttp,fpi,dntp)
    
if __name__=='__main__':
    LV=[]   
    TPS=[]
    DNTP=[]
    auconly=False # whether to calculate the avg. auc of the complexes or do more
    (auc,(fp,tp),(A,Rx,D,L,cids,r,dkey))=getAUC('result_n.res.pkl')
    dAo=dict(zip(cids,A)) #AUCs from the training data set (CV) only
    loopath='./DBD3LOOCV/'
    fs=glob.glob(loopath+'*.pairpred.txt')
    dA={}
    dsA={}
    for i,ifile in enumerate(fs):
        cid=getFileParts(getFileParts(ifile)[1])[1]
        print 'cid =',cid,100*float(i+1)/len(fs),'% done'
        #
        if auconly:
            auc=readFile(ifile,auconly=True)
        else:
            (auc,Mv,Ml,lseq,rseq,lrV,rrV)=readFile(ifile);#(auc,Mv,Ml,lseq,rseq,lrV,rrV)

            #pdb.set_trace()
            (la,lv,ll)=getAUC4Protein(lrV)
            (ra,rv,rl)=getAUC4Protein(rrV)
            Mvx=Mv.ravel()
            Mlx=Ml.ravel()
            nidx=~np.isnan(Mvx) &  ~np.isnan(Mlx)
            Mvx[~nidx]=-np.inf            
            (ttp,fpi,dntp)=findNTPinTop(Mvx,Mlx,Mv.shape,top=200)
            Mvx=Mvx[nidx]
            Mlx=Mlx[nidx]
            pp=np.sum(Mlx==1) # total number of positives
            nn=len(Mlx)-pp #total number of negatives
            TPS.append([ttp,fpi,100.0*ttp/pp,pp,nn,pp+nn])
            DNTP.append(dntp)
            #LV.append((list(Mvx),list(Mlx)))       
            dsA[cid]=(auc,la,ra)
        #pdb.set_trace()
        dA[cid]=auc
        print cid,auc,dAo[cid]
    print 'Complex wise AUC = ',np.mean(dA.values()),'AUC for reduced set = ',np.mean([dAo[k] for k in dA.keys()])
    if not auconly:
        p12=map(list,zip(*dsA.values()));pa=p12[0];p1=p12[1];p2=p12[2];ps=p1;ps.extend(p2);
        print 'Complex Wise AUC =',np.mean(pa),'Protein Wise AUC =',np.mean(ps)
    plt.hist(np.array(DNTP).flatten(),[0,1,2,3,4,5,6,1000],cumulative=True);plt.grid();plt.xlabel('sequence distance');plt.ylabel('counts');plt.title('Number of top 200 predictions vs. sequence distance from nearest true positive');plt.show()
    [np.sum(dn<2.0) for dn in DNTP]
    cids=[getFileParts(getFileParts(ifile)[1])[1] for ifile in fs]
    [dsA[cid] for cid in cids]
    [dAo[cid] for cid in cids]
    DD=[]
    for dn in DNTP:
        d=np.nonzero(dn<2.0)[0]
        if len(d):
            DD.append(d)
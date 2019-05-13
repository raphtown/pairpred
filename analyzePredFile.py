# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 19:50:07 2013
Code for displaying and analyzing the results from PAIRPred
@author: root
Usage:
    python analyzePredFile.py -i ./Flaviviridae/NS3_NS5_DBD4.InterPred.txt -o tmp.txt -n 10 -w 5 -l A -r A -t 2

"""
import numpy as np
#import matplotlib
#matplotlib.use('wx')
import matplotlib.pyplot as plt
from copy import deepcopy
import re
import sys
from optparse import OptionParser

def parseLine(lns):
    lpid=lns[0].split(':')
    if len(lpid)<2:
        lcid=''
    else:
        lcid=lpid[1]    
    lid=int(re.sub("\D", "", lns[2]))
    rpid=lns[5].split(':')
    if len(rpid)<2:
        rcid=''
    else:
        rcid=rpid[1]    
    rid=int(re.sub("\D", "", lns[7])) #re.sub("\D", "", "aas30dsa20")
    lrid=int(lns[4])
    rrid=int(lns[9])
    v=float(lns[10])
    l=int(lns[11])
    ls=lns[1]
    rs=lns[6]
    
    return (lcid,lid,rcid,rid,v,l,ls,rs,lrid,rrid)
    
def readFile(ifile,ilcid=None,ircid=None,auconly=False):
    """
    Read result file (ifile) produced by PAIRPred and return (auc,Mo,lseq,rseq)
        auc: AUC score -- None if it does not exist in the file
        Mo: Matrix of scores with the rows corresponding to the ligand and 
            the columns corresponding to the receptor. The index of an element
            indicates the location of the residues in the PDB file
        lseq: ligand sequence
        rseq: Receptor sequence
    """
    auc=None
    Lid=[]
    Rid=[]
    V=[]
    L=[]
    lseq={}
    rseq={}
    lrV={}
    rrV={}
    for ln in open(ifile,'r'):
        lns=ln.split()
        if lns[0]=='#': #ignore comments
            continue
        elif lns[0]=='AUC':
            auc=float(lns[2])
            if auconly:
                return auc
        else:
            (lcid,lid,rcid,rid,v,l,ls,rs,lrid,rrid)=parseLine(lns)
            if lrid not in lrV:
                lrV[lrid]=(v,l)                
            else:
                lrV[lrid]=(np.max((lrV[lrid][0],v)),np.max((lrV[lrid][1],l)))
            if rrid not in rrV:
                rrV[rrid]=(v,l)                
            else:
                rrV[rrid]=(np.max((rrV[rrid][0],v)),np.max((rrV[rrid][1],l)))            
      
            if (ilcid is None and ircid is None) or (lcid==ilcid and rcid==ircid):
                Lid.append(lid)
                Rid.append(rid)
                V.append(v)
                L.append(l)
                lseq[lid]=ls
                rseq[rid]=rs
    #
    for ri in range(np.max(Lid)):
        if ri not in lseq:
            lseq[ri]='-'
    for ci in range(np.max(Rid)):
        if ci not in rseq:
            rseq[ci]='-'    
    assert len(np.unique(np.diff(lseq.keys())))==1 and (np.unique(np.diff(lseq.keys())))[0]==1
    assert len(np.unique(np.diff(rseq.keys())))==1 and (np.unique(np.diff(rseq.keys())))[0]==1
    Mv=np.zeros((np.max(Lid)+1,np.max(Rid)+1))
    Mv.fill(np.nan)
    Ml=np.zeros(Mv.shape)
    Ml.fill(np.nan)
    #pdb.set_trace()
    lseq=lseq.values();lseq=''.join(lseq)
    rseq=rseq.values();rseq=''.join(rseq)     
    for i in range(len(Lid)):
        Mv[Lid[i],Rid[i]]=np.nanmax((Mv[Lid[i],Rid[i]],V[i]))
        Ml[Lid[i],Rid[i]]=np.nanmax((Ml[Lid[i]-1,Rid[i]],L[i]))
    return (auc,Mv,Ml,lseq,rseq,lrV,rrV)
    
def plotPred(Mo,thr=2.5,lname='L',lseq='',rname='R',rseq='',met='PAIRPred'):
    """
    Plots the results from PAIRPred. It displays the matrix of scores as an image
    It also generates a imahe in which only the top 'thr' % of scores are displayed.
    It generates the plots for scores for the ligand and receptor proteins as well.
        Mo: score matrix
        thr: select top thr % only
        lname: name of the ligand
        lseq: sequence of the ligans
        rname: ma,e pf tje receptor
        rseq: Seqiemce of the receptor
        met: Name of the method
    """
    xstep=10
    ystep=10
    met=met+': '
    cmap =plt.cm.jet;# plt.get_cmap()
    cmap.set_bad(color = 'k', alpha = 1.)
    Mf=Mo.flatten()
    Mf=Mf[~np.isnan(Mf)]
    plt.figure()
    (c,b,p)=plt.hist(Mf,bins=1000,cumulative=True,normed=True);#plt.xlabel('value');plt.ylabel('Normalized Frequency');plt.title(met+'distribution')
    plt.close()
    thrv=b[np.where(c<(1-thr/100.0))[0][-1]]
    #xtn=range(0,Mo.shape[1],50);xt=[Lpid[x] for x in xtn]
    plt.figure();plt.imshow(Mo,cmap=cmap);plt.colorbar();plt.xlabel(rname);plt.ylabel(lname);plt.title(met+'All Predictions');plt.xticks(range(0,Mo.shape[1],xstep));plt.yticks(range(0,Mo.shape[0],ystep));
    M1=deepcopy(Mo)
    M1[Mo<thrv]=np.nan
    
    plt.figure();plt.imshow((M1),cmap=cmap);plt.colorbar();plt.xlabel(rname);plt.ylabel(lname);plt.title(met+'Top '+str(thr)+'% Predictions');plt.xticks(range(0,Mo.shape[1],xstep));plt.yticks(range(0,Mo.shape[0],ystep));
    plt.figure();plt.plot(np.nanmax(Mo,axis=0));plt.title(met+rname+' Prediction');plt.xticks(range(0,Mo.shape[1],xstep));plt.grid()
    plt.figure();plt.plot(np.nanmax(Mo,axis=1));plt.title(met+lname+' Prediction');plt.xticks(range(0,Mo.shape[0],ystep));plt.grid()
def sortScores(Mo):
    M1=deepcopy(Mo)
    Mor=M1.ravel()
    nidx=np.isnan(Mor)
    Mor[nidx]=-np.inf
    Mori=np.argsort(Mor)[::-1]
    Mori=[i for i in Mori if Mor[i]!=-np.inf]
    (r,c)=np.unravel_index(Mori,Mo.shape)
    v=Mor[Mori]
    
    return (r,c,v)
def writeTop2File(Mo,fname,lseq='',rseq='',HWS=5,top=np.inf):
    """
    Given the score matrix, it writes the top 'top' pairs into a file
    178<\t>	247<\t>1.43149930273<\t>EPDYEVDEDIF<\t>	V<\t>	RFTMRHKKATY<\t>H
    ligand_resid,receptor_resid,score,HWS sized sequnce window around the residue, residue...
    """
    
    def getslseq(lseq,HWS,ri):
        rle=ri-HWS
        rre=ri+HWS+1
        rlpad=rrpad=''
        if rle<1:                
            rlpad='*'*(0-rle+1)
           #pdb.set_trace()
            rle=1            
        if rre>len(lseq)-1:
            rrpad='*'*(rre-len(lseq)+1)
            #pdb.set_trace()
            rre=len(lseq)-1            
        slseq=rlpad+lseq[rle:rre]+rrpad
        return slseq    
        
    (r,c,v)=sortScores(Mo)
    f = open(fname, 'w+')
    s='#lig_id\trec_id\tscore\tlig_res\tlig_win\trec_res\trec_win\n'
    f.write(s)
    for i in range(int(np.min((len(v),top)))):
        s=str(r[i]) + '\t'+str(c[i])+'\t'+str(v[i])
        if len(lseq):
            slseq=getslseq(lseq,HWS,r[i])
            s=s+'\t'+lseq[r[i]]+'\t'+slseq
        else:
            s=s+'\t'+'-'
        s=s+'\t'
        if len(rseq):
            srseq=getslseq(rseq,HWS,c[i])
            s=s+'\t'+rseq[c[i]]+'\t'+srseq
        else:
            s=s+'\t'+'-'      
        
        f.write(s+'\n')
    f.close()
    #return (r,c,v)
if __name__=="__main__":
    #ifile,ilcid,ircid,topfname,top,HWS,plotthr
    parser = OptionParser()
    parser.add_option("-i", "--ifile", dest="ifile",
                      help="Input Prediction file", metavar="file path")
    parser.add_option("-l", "--lchain", dest="ilcid",
                      help="Ligand chain (default: All)", metavar="character")
    parser.add_option("-L", "--lname", dest="lname",
                      help="Ligand name (default: L)", metavar="character")         
    parser.add_option("-R", "--rname", dest="rname",
                      help="Receptor name (default: R)", metavar="character")                          
    parser.add_option("-r", "--rchain", dest="ircid",
                      help="Receptor chain (default: All)", metavar="character")
    parser.add_option("-o", "--ofile", dest="ofile",
                      help="Output File", metavar="file path")
    parser.add_option("-n", "--num", dest="num",
                      help="number of top predictions in output file (default: 100)", metavar="Number")    
    parser.add_option("-w", "--win", dest="win",
                      help="sequence half-window size (default: 5)", metavar="Number")             
    parser.add_option("-t", "--thr", dest="thr",
                      help="top % threshold in plotting (default: 1 %)", metavar="Percentage number")                         
    (options, args) = parser.parse_args()
    
    if options.ifile is None or options.ofile is None:
        print 'Must specify both input and output file names. Exiting.'
        print 'Type python analyzePredFile.py -h for help'
        sys.exit(1)
    if options.lname is None:
        options.lname='L'
    if options.rname is None:
        options.rname='R'        
    if options.win is None:
        options.win=5
    else:
        options.win=int(options.win)
    if options.thr is None:
        options.thr=1.0
    else:
        options.thr=float(options.thr)
    if options.num is None:
        options.num=100
    else:
        options.num=int(options.num)    
    print options
    (auc,Mv,Ml,lseq,rseq,lrV,rrV)=readFile(options.ifile,ilcid=options.ilcid,ircid=options.ircid)
    if auc is not None:
        print 'AUC =', auc
    writeTop2File(Mv,fname=options.ofile,lseq=lseq,rseq=rseq,HWS=options.win,top=options.num)
    plotPred(Mv,lname=options.lname,rname=options.rname,met='PAIRPred',thr=options.thr)
    plt.show()
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 13:32:18 2013
Script for generating and visualizing the links (.lnx) files from a PAIRPred
or other predictor
@author: Afsar
"""
import os
from Bio.PDB import *
from analyzePredFile import readFile,sortScores

def getFileParts(fname): #Added to decouple it from utils
    "Returns the parts of a file"
    (path, name) = os.path.split(fname)
    n=os.path.splitext(name)[0]
    ext=os.path.splitext(name)[1]
    return (path,n,ext)
def readLinkFile(fname,top=0): 
    """
    Reader for lnx files
    top=0 implies that all the file is read. Otherwise only the top 'top'
    links are read
    """
    L=[]
    lpdb=None
    lmol=None
    rpdb=None
    rmol=None
    cnt=0
    for ln in open(fname,'r'):
        lns=ln.split()
        if lns[0]=='#':
            continue
        if lns[0]=='lpdb':
            lpdb=lns[1]
            (_,lmol,_)=getFileParts(lpdb)
            continue
        if lns[0]=='rpdb':
            rpdb=lns[1]
            (_,rmol,_)=getFileParts(rpdb)
            continue
        cnt=cnt+1
        L.append(((lns[0],int(lns[1])),(lns[2],int(lns[3])),float(lns[4])))
        if top and top==cnt:
            break
    return ((lpdb,lmol),(rpdb,rmol),L)
    
def writeLinkFile(res,lpdb,rpdb,ldmap,rdmap,ofname=None):
    """
    Writes the links files.
    res: It can either be a result file or an interaction matrix
       result file name, e.g. ./BDB3LOOCV/1N2C.pairpred.txt
    lpdb: ligand pdb file (bound state)
    rpdb: receptor pdb file (bound state)
    lpdb and rpd can be same if the complex is containted in a single file
    These files can be different from the ones used to generate 'res'
    So we require the mapping between the chains in res and in the lpdb,rpdb
    ldmap: Mapping between the molecules and chains in res to those in lpdb
        e.g. ldmap=('A','C') means that chain A in res corresponds to chain C in lpdb
    rdmap: Same as ldmap except that it is for the receptor
    ofname: output file name [Default: 'set'+ldmap[0]+rdmap[0]+'_'+ldmap[1]+rdmap[1]+'.lnx']    
    """    
    res_lc,pdb_lc=ldmap   
    res_rc,pdb_rc=rdmap
    if ofname is None:
        ofname='set'+res_lc+res_rc+'_'+pdb_lc+pdb_rc+'.lnx'
    if type(res)==type(''):
        (_,Mv,_,_,_,_,_)=readFile(res,ilcid=res_lc,ircid=res_rc)
    else:
        Mv=res
    (r,c,v)=sortScores(Mv)
    f=open(ofname,'w+')
    s='lpdb '+lpdb;    f.write(s+'\n')
    s='rpdb '+rpdb;    f.write(s+'\n')
    for i in range(len(v)):
        s=pdb_lc+' '+str(r[i])+' '+pdb_rc+' '+str(c[i])+' '+str(v[i]);    f.write(s+'\n')
    f.close()
def launchPymol(lpdb,rpdb,flist=[],top=100):
    """
    Launches Pymol and initiates the commands addLinks and clrLinks
    lpdb: ligand pdb file
    rpdb: receptor pdb file
    flist: [optional] is a list of .lnx files that need to be displayed
    top,dthr are used only when flist is specified.
    """    
    import __main__    
    __main__.pymol_argv = [ 'pymol']
    import pymol
    from pymol import cmd
    pymol.finish_launching()           
    
    cmd.load(lpdb)
    if lpdb!=rpdb:
        cmd.load(rpdb)
    cmd.show('cartoon',"all")
    cmd.hide('lines',"all")
    cmd.do('util.cbc')
    cmd.do('run showLinksInPymol.py,main ')        
    for f in flist:
        cmd.do('loadLink '+f+',' +str(top))
    #cmd.do('showLinks '+str(dthr))
if __name__=='__main__':
    
    
    lpdb='./testcase/DBD4_IGS15_NS1B_PairPred/3SDL.pdb'    
    rpdb='./testcase/DBD4_IGS15_NS1B_PairPred/3SDL.pdb'  
    """
    resfname='./testcase/DBD4_IGS15_NS1B_PairPred/HNS1_dbd4.InterPRed.txt'
    flist=['./setCA_PAIRpred.lnx' ,'./setDB_PAIRpred.lnx','./setCB_PAIRpred.lnx','./setDA_PAIRpred.lnx']
    ldmap=('A','C')
    rdmap=('B','A')
    lnkfile=flist[0]
    writeLinkFile(resfname,lpdb,rpdb,ldmap,rdmap,ofname=lnkfile)     
    ldmap=('A','D')
    rdmap=('B','B')
    lnkfile=flist[1]
    writeLinkFile(resfname,lpdb,rpdb,ldmap,rdmap,ofname=lnkfile)
    ldmap=('A','C')
    rdmap=('B','B')
    lnkfile=flist[2]
    writeLinkFile(resfname,lpdb,rpdb,ldmap,rdmap,ofname=lnkfile)
    ldmap=('A','D')
    rdmap=('B','A')
    lnkfile=flist[3]
    writeLinkFile(resfname,lpdb,rpdb,ldmap,rdmap,ofname=lnkfile)  
    """
    flist=['./testcase/DBD4_IGS15_NS1B_PairPred/setCA_PAIRpred200.lnx',\
    './testcase/DBD4_IGS15_NS1B_PairPred/setCB_PAIRpred200.lnx',\
    './testcase/DBD4_IGS15_NS1B_PairPred/setDA_PAIRpred200.lnx',\
    './testcase/DBD4_IGS15_NS1B_PairPred/setDB_PAIRpred200.lnx']
    launchPymol(lpdb,rpdb,flist,200)
    
    """
    resfname='./DBD3LOOCV/2UUY.pairpred.txt'
    lpdb='./DBD4N/DBD4/2UUY_l_b.pdb'    
    rpdb='./DBD4N/DBD4/2UUY_r_b.pdb'    
    ldmap={'A':'B'}
    rdmap={'A':'A'}
    lnkfile='./tmp.lnx'
    writeLinkFile(resfname,lpdb,rpdb,ldmap,rdmap,ofname=lnkfile) 
    launchPymol(lpdb,rpdb)
    """
    
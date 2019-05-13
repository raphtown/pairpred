# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 16:31:29 2012
@author: Afsar
This module extracts information from the bound structures involved in a complex
It performs mapping between the sequences of the bound and unbound sequences
which are extracted from their structure files. This is done using global
sequence alignment. This mapping is then used in extracting positive and
negative training examples in getExamplesDBD.py
"""
from BISEPutils import *
from myPDB import *
import Bio.pairwise2

def mapU2B(us,uS2Ri,ulR,bs,bS2Ri,blR):
    """
    Get the mapping of indices between the bound and unbound residues
    us: sequence of the unbound protein
    uS2Ri: seq. to R index list for unbound protein
    ulR: length of R list for the unbound protein
    Similar for the bound protein
    Return: u2b,b2u
    #uR[i] corresponds to bR[u2b[i]] if u2b[i] is not nan
    #bR[i] corresponds to uR[b2u[i]] if b2u[i] is not nan
    """
    u2b=[np.nan for k in range(ulR)]
    b2u=[np.nan for k in range(blR)]
    aln=Bio.pairwise2.align.globalxs(us,bs,-1,-0.1)[0]
#    print aln
    ui=0
    bi=0
    for k in range(len(aln[0])):
        uc=aln[0][k]
        bc=aln[1][k]
        if uc==bc:
            u2b[uS2Ri[ui]]=bS2Ri[bi]
            b2u[bS2Ri[bi]]=uS2Ri[ui]
        ui=ui+(uc!='-')
        bi=bi+(bc!='-')
    pu2b=sum(np.isnan(u2b))/float(len(u2b))
    pb2u=sum(np.isnan(b2u))/float(len(b2u))
#    if  pu2b > 0.01:
#        print "Warning: "+str(pu2b*100)+" % of unbound residues not matched"
#    if  pb2u > 0.01:
#        print "Warning: "+str(pb2u*100)+" % of bound residues not matched"
    return (u2b,b2u)

class myPDBComplex:
    """
    Class responsible for extracting information from the bound structures of 
    the proteins.
    Attributes:
        fname: list of paths of files constituting the complex
        N: length of fnames
        cid: list of list of chain ids of each file
        R: R[i] is the list of residues (biopython) for file i
        Coords: Coords[i][j] is a list of coordinates of residue R[j] for file i        
        seq: seq[i] contains the combination of the peptide sequence
        dthr: distance threshold
        D: D[i][j] is the list of distances between files i and j
    """
    def __init__(self,fnames,dthr=6.0):
        
        self.fnames=fnames
        self.N=len(fnames)
        self.dthr=dthr
        
        self.Coords=[]
        self.R=[]
        self.D=[[[] for i in range(self.N)] for j in range(self.N)]
        self.seq=[]
        #self.stx=[]
        self.S2Ri=[]
        for (i,f) in enumerate(fnames):
            (_,R,_,seq,S2Ri)=readPDB(f)      #(cid,stx,R,pp,seq,S2Ri)      
#            self.cid.append(cid)
            #self.stx.append(stx)
            self.seq.append(seq)
            self.R.append(R)
            self.S2Ri.append(S2Ri)
            self.Coords.append(getCoords(R))            
            for j in range(0,i):
                (self.D[i][j],self.D[j][i])=getDist(self.Coords[i],self.Coords[j],dthr)
        
        """
        #Required for ASA computations
        out_file = tempfile.NamedTemporaryFile(suffix='.pair.pdb')
        out_file.close()
        tmpfile=out_file.name
        mergePDBFiles(fnames,tmpfile)
        
        os.remove(tmpfile)
        """
    def findNforU(self,ufnames):
        """
        returns the distance object but changes the indices to represent the 
        mapping between bound and unbound structures
        """
        D=[[[] for i in range(self.N)] for j in range(self.N)]
        U2B=[]
        B2U=[]
        LenU=[]
        for (idx,f) in enumerate(ufnames):
            (_,uR,_,us,uS2Ri)=readPDB(f)
            (bR,bs,bS2Ri)=(self.R[idx],self.seq[idx],self.S2Ri[idx])
            (u2b,b2u)=mapU2B(us,uS2Ri,len(uR),bs,bS2Ri,len(bR))
            LenU.append(len(uR))
            U2B.append(u2b)
            B2U.append(b2u)
        for i in range(self.N):
            for j in range(i+1,self.N):
                oD=self.D[i][j]
                for k in range(len(oD)):
                    ov=oD[k]
                    vij=(B2U[i][ov[0]],B2U[j][ov[1]],ov[2])
                    vji=(B2U[j][ov[1]],B2U[i][ov[0]],ov[2])
                    D[i][j].append(vij)
                    D[j][i].append(vji)
        return (D,LenU)
if __name__=="__main__":
    #SIMPLE TESTING CODE
    pdbpath='.//benchmark4//structures//Benchmark_4_updated//'    
    pdbid="1GL1"
    lbfname=pdbpath+pdbid+"_"+"l_b.pdb"
    rbfname=pdbpath+pdbid+"_"+"r_b.pdb"
    lufname=pdbpath+pdbid+"_"+"l_u.pdb"
    rufname=pdbpath+pdbid+"_"+"r_u.pdb"
    
    C=myPDBComplex([lbfname,rbfname])
    (D,_)=C.findNforU([lufname,rufname])
    (_,ruR,_,_,_)=readPDB(rufname)
    (_,luR,_,_,_)=readPDB(lufname)
    print len(D[0][1])
    def test(idx):       
        x=(C.D[0][1][idx],D[0][1][idx]);
        print x;
        print getDist(getCoords([C.R[0][x[0][0]]]),getCoords([C.R[1][x[0][1]]]),thr=np.inf);
        print getDist(getCoords([luR[x[1][0]]]),getCoords([ruR[x[1][1]]]),thr=np.inf);
        print C.R[0][x[0][0]],luR[x[1][0]],C.R[1][x[0][1]],ruR[x[1][1]]
        print C.R[0][x[0][0]].get_full_id(),luR[x[1][0]].get_full_id(),C.R[1][x[0][1]].get_full_id(),ruR[x[1][1]].get_full_id()

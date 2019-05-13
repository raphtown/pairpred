# -*- coding: utf-8 -*-
"""
Created on Thu Oct 04 01:30:15 2012
@author: Afsar
This module contains helper functions.
"""
msmspath='/scratch/PI/rondror/ppi/pairpred/lib/msms' #contains stride and msms implementation
stridepath='/scratch/PI/rondror/ppi/pairpred/lib/stride'
blastpath='/scratch/PI/rondror/ppi/pairpred/lib/psiblast/ncbi-blast-2.6.0+/bin'
blastdbpath='/scratch/PI/rondror/ppi/non-redundant'
spxpath='/scratch/PI/rondror/ppi/pairpred/lib/spinex/spineXpublic/' #spinex path
psaiapath='/scratch/groups/rondror/ppi/pairpred/lib/psaia-src/PSAIA_1.0_source/bin/linux/psa' #psaia path

import sys
import matplotlib
matplotlib.use('Agg') #Must be uncommented for running on department machines
import os
if psaiapath not in  os.environ["PATH"]:
    print 'Adding '+psaiapath+' to system path'
    os.environ["PATH"]+=os.pathsep+psaiapath
if msmspath not in  os.environ["PATH"]:
    print 'Adding '+msmspath+' to system path'
    os.environ["PATH"]+=os.pathsep+msmspath
if stridepath not in  os.environ["PATH"]:
    print 'Adding '+stridepath+' to system path'
    os.environ["PATH"]+=os.pathsep+stridepath
if blastpath not in  os.environ["PATH"]:
    print 'Adding '+blastpath+' to system path'
    os.environ["PATH"]+=os.pathsep+blastpath
if spxpath not in  os.environ["PATH"]:
    print 'Adding '+spxpath+' to system path'
    os.environ["PATH"]+=os.pathsep+spxpath
import pdb
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
import cPickle
from Bio.PDB import *
from Bio.PDB.Polypeptide import *
from Bio import SeqIO
import tempfile
from scipy.sparse import lil_matrix
from Bio.SubsMat import MatrixInfo

AA='ACDEFGHIKLMNPQRSTVWY-'
aaidx=dict(zip(AA,range(len(AA))))

def score_match(pair,matrix=MatrixInfo.blosum62):
    """
    Given a tuple pair of amino acids, it returns the substitution matrix
    score
    """

    if pair not in matrix:
        pair=tuple(reversed(pair))
    if pair not in matrix:
        return 0.0
    else:
        return matrix[pair]
def getSubMat(matrix=MatrixInfo.blosum62):
    """
    Returns a dictionary representation of the columns of a substitution matrix
    """
    M={}
    for a in AA:
        M[a]=[]
        for b in AA:
            M[a].extend([score_match((a,b),matrix)])
    return M

BLOSUM62=getSubMat()

def getWPSSM(xpm,HWS=5):
    """
    Given a np.array xpm this function creates the windowed representation
    (for use in PSSM, PSFM etc)
    """
    HWS=int(HWS)
    (d,N)=xpm.shape
    pm=np.hstack((np.zeros((d,HWS)),xpm,np.zeros((d,HWS))))
    ws=2*HWS+1
    wpm=np.zeros((ws*d,N))
    for i in range(N):
        wpm[:,i]=pm[:,i:i+ws].flatten('F')
    return wpm

def getSubMatFeats(s,HWS=5):
    """
    Given a sequence, this function returns the subsitution matrix representation
    """
    HWS=int(HWS)
    smat=np.zeros((len(AA),len(s)))
    for (i,a) in enumerate(s):
        try:
            smat[:,i]=BLOSUM62[a]
        except Exception as e:
            print e
            continue
    return getWPSSM(smat,HWS)


def getDSSP(stx,fname):
    """
    Biopython's dssp does not process broken chains (or missing residues). So
    what we do here is to write a temp pdb file for each peptide and apply DSSP
    on it. The return is a dictionary object in which the DSSP proeprties for
    all residues for all chains have been merged.
    """
    class pepSelect(Select):
        """
        Required for selecting the residues within a peptide to write them to a
        file.
        """
        def __init__(self,p):
            self.pL=[r.get_full_id() for r in p]
        def accept_residue(self,res):
            #pdb.set_trace()
            resid=res.get_full_id()#[:-1]
            if resid  in self.pL:
                return 1
            else:
                return 0
    pp=PPBuilder().build_peptides(stx[0])
    io=PDBIO()
    io.set_structure(stx)
    dssp=dict()
    out_file = tempfile.NamedTemporaryFile(suffix='.dssp')
    tmpfname=out_file.name
    out_file.close()
    try:
        for p in pp:
            io.save(tmpfname,pepSelect(p))
            d=DSSP(stx[0], tmpfname)
            dssp=dict(dssp,**d)
    except:
            e = sys.exc_info()[0]
            print e
            pdb.set_trace()
            print  "Oops! Problem running DSSP! Is it installed correctly?"
    finally:
        os.remove(tmpfname)
    return dssp

def getFileParts(fname):
    "Returns the parts of a file"
    (path, name) = os.path.split(fname)
    n=os.path.splitext(name)[0]
    ext=os.path.splitext(name)[1]
    return (path,n,ext)
def getResiId(fid):
    """
    Given the full id of a residue, return the tuple id form
    """
    (_,_,cid,(_,ridx,rinum))=fid
    return (cid,str(ridx)+rinum.strip())#

def getResLetter(r2):
    """
    Get the letter code for a biopython residue object
    """
    r2name=r2.get_resname()
    if is_aa(r2name, standard=True):
        scode=three_to_one(r2name)
    else:
        scode='-'
    return scode

def getSideChainV(r):
    """
    Find the average of the unit vectors to different atoms in the side chain
    from the c-alpha atom. For glycine the average of the N-Ca and C-Ca is
    used.
    Returns (C-alpha coordinate vector, side chain unit vector) for residue r
    """
    u=None
    gly=0
    if is_aa(r) and r.has_id('CA'):
        ca=r['CA'].get_coord()
        dv=np.array([ak.get_coord() for ak in r.get_unpacked_list()[4:]])
        if len(dv)<1:
            if r.has_id('N') and r.has_id('C'):
                dv=[]
                dv.append(r['C'].get_coord())
                dv.append(r['N'].get_coord())
                dv=np.array(dv)
                gly=1
            else:
                #pdb.set_trace()
                return None
        dv=dv-ca
        if gly:
            dv=-dv
        n=np.sum(np.abs(dv)**2,axis=-1)**(1./2)
        v=dv/n[:,np.newaxis]
        u=(Vector(ca),Vector(v.mean(axis=0)))
    return u

def getCoords(R):
    """
    Get atom coordinates given a list of biopython residues
    """
    Coords=[]
    for (idx,r) in enumerate(R):
        v=[ak.get_coord() for ak in r.get_list()]
        Coords.append(v)
    return Coords

def getDist(C0,C1,thr=np.inf):
    """

    """
    N0=[]
    N1=[]
    for i in range(len(C0)):
        for j in range(len(C1)):
            d=spatial.distance.cdist(C0[i], C1[j]).min()
            #dji=spatial.distance.cdist(C1[j], C0[i]).min()
            #d=min(dij,dji)
            #print d
            if(d<thr):# and not np.isnan(self.Phi[i]) and not np.isnan(self.Phi[j])
                N0.append((i,j,d))
                N1.append((j,i,d))
    return (N0,N1)
def readFASTA(fname):
    """
    Reads the fasta file fname and returns the sequence only
    """
    handle = open(fname, "rU")
    record = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    if len(record)>1:
        print "Warning: Input FASTA file must have only one protein sequence in it. Using only the first sequence."
    #pdb.set_trace()
    record=record[0]
    return record.seq.tostring().upper()

def readPDB(fname,name=None):
    """
    Extract info from a PDB file
        fname: path of pdb file
        name: name of the structure (default name of the file without extension)
        return:: (stx,R,pp,seq,S2Ri)

            stx: structure object
            R: list of residues
            pp: list of polypeptides in the structure
            seq: combined sequence (for all polypeptides)
            S2Ri: Sequence to R mapping index list, seq[i] corresponds to
                R[S2Ri[i]]
    """
    if name is None:
        (_,name,_)=getFileParts(fname)
    stx=PDBParser() .get_structure(name,fname)
    assert len(stx)==1 #there should be only one structure
#    cid=[]
#    for c in stx[0].get_list():
#        cid.append(c.id)
    R=Selection.unfold_entities(stx,'R') #list of residues
    pp=PPBuilder().build_peptides(stx)
    seq=''.join([p.get_sequence().tostring() for p in pp])
    Rdict=dict(zip(R,range(len(R))))
    S2Ri=[Rdict[r] for p in pp for r in p]
    return (stx,R,pp,seq,S2Ri)

def getSeqFV(stx,R,HWS=10):
    """
    Return a FV based on sequence (uses structure)
    """
    HWS=int(HWS)
    #FV=[[] for x in range(len(R))]
    FV=None
    Rdict=dict(zip(R,range(len(R))))
    first=1
    #pdb.set_trace()
    for c in stx[0]:
        pp=PPBuilder().build_peptides(c)
        if len(pp)==0:
            print "Ignored the empty chain encounted in ",stx.get_full_id()
            continue
        (s,s2r)=getSeqForChain(pp,Rdict)
        s='-'*HWS+s+'-'*HWS
        idx=0
        for m in range(HWS,len(s)-HWS):
            if s[m]!='-':
                w=s[m-HWS:m+HWS+1]
                (pd1,p1)=getPD1Spec(w)
                #pdb.set_trace()
                fv=np.vstack([np.array(pd1.todense()),np.array(p1.todense())])
                if(first):
                    first=0
                    #create k based on legnth of the vector
                    FV=lil_matrix((fv.shape[0],len(R)))
                FV[:,s2r[idx]]=fv
                idx=idx+1
    FVF=FV.todense()
    #If no features for the whole residue, set to nan
    FVF[:,np.nonzero(np.array(np.sum(FVF,0)==0).ravel())[0]]=np.nan
    FV=lil_matrix(FVF)
    #pdb.set_trace()
    return FV
def getSeqForChain(pp,Rdict):
    """
    For a single chain only!!!
    Get the sequence for a chain and the S2Ri associated with non '-' chars in s
    Remember len(s) is equal to len(S2Ri) iff there are no dashes in s
    """
    s=[]
    s.append(pp[0].get_sequence().tostring())
    for idx in range(len(pp)-1):
        b=pp[idx][-1].id[1]+1
        e=pp[idx+1][0].id[1]-1
        d='-'*(e-b+1)
        s.append(d)
        s.append(pp[idx+1].get_sequence().tostring())
    S2Ri=[Rdict[r] for p in pp for r in p]
    s=''.join(s)
    return (s,S2Ri)

def getPD1Spec(s,param=None):
    """
    Get 1-spectrum representation of s (ignoring '-'')
    """
    dv=np.sqrt(1/20.0)
    V = np.zeros((len(AA),len(s)),dtype='float64')
    for k in range(len(s)):
        if s[k]!='-':
            try:
                V[aaidx[s[k]],k]=V[aaidx[s[k]],k]+1.0
            except KeyError:
                pass
        else:
            for a in aaidx:
                if a!='-':
                    V[aaidx[a],k]=V[aaidx[a],k]+dv
    #pdb.set_trace()
    #m=np.sqrt(len(s))
    v=V.reshape((np.prod(V.shape),1))
    v=lil_matrix(v/np.linalg.norm(v))
    v1=V.sum(axis=1)
    v1=lil_matrix(v1/np.linalg.norm(v1)).T
    return (v,v1)

def mergePDBFiles(fnames,ofname):
    ofh=open(ofname,'w')
    for f in fnames:
        fh=open(f,'r')
        d=fh.read()
        fh.close()
        ofh.write(d)
    ofh.close()

def copy_dict(d, *keys):
    """Make a copy of only the `keys` from dictionary `d`."""
    return {key: d[key] for key in keys}

def chunks(l, n):
    """ Yield successive n-sized chunks from list l, returns list of lists."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]
def mergeDicts(Dlist):
    """ merges a list of dictionaries into a single dictionary"""
    dld={}
    for d in Dlist: dld.update(d)
    return dld
def combineList(l):
    return [item for sublist in l for item in sublist]

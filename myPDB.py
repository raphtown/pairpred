# -*- coding: utf-8 -*-
"""
Created on Tue Oct 02 21:06:18 2012

@author: Afsar
"""

from BISEPutils import *
from scipy.spatial.distance import *
from scipy import linalg as LA
from psiblast import *
import stride
import dill
import spinex
import psaia
import os
import glob
class myPDB:
    """
    Version: 1240HRS30DEC12
    Author: Fayyaz
    Given a PDB file, this class is responsible for extracting residue level
    information to be used by 'mySimMtx' to evaluate the kernels
    Data members:
        name: name of the structure (default: file name of the PDB file)
        ifname: input file path
        stx: biopython structure object
        R: List of residues (from biopython) in the filed
        pp: list of polypeptides
        seq: combined sequnce
        S2Ri: sequrnce to R index, seq[i] correspoinds to residue R[S2R[i]]
        Coords: Type list with each element of the list representing a residue
            and containing a list of atomic coordinates of that residue
        S: Self-Similarity matrix (see getSimMtx)
        sg: sigma value used in S
        thr: threhsold used in S
        UC: aa composition in the up direction of a residue (inc. res. itself)
        UN: aa count in the up direction of a residue
        DC: aa composition in the down direction of a residue (inc. res. itself)
        DN: aa count in the down direction of a residue
        CN: Coordination number of a residue ~(UN+DN)
        B: Numpy array with B[ridx] containing the maximum B-factor of any atom
            in the residue R[ridx]
        SS: Secondary structure representation for each residue
        ASA: Accessible surface area for each residue
        rASA: relative accessible surface area for each residue
        asa: ASA predicted from sequence using spinex (requires auxdir)
        rasa: rasa predicted from sequence
        pssm: PSSM returned by BLASTPGP
        psfm: PSFM returned by BLASTPGP
        info: Information content returnd by BLASTPGP
        wpssm: Windowed pssm
        wpsfm: Windowed psfm
        wbsmf: Windowed Blosum62 substitution matrix features
        HWS_seqFV: HWS for sequence features
        HWS_PSSM: HWS for wpssm,wpsfm,wbsmf
        Phi: Phi torsion angle for each residue
        Psi: Psi torsion anfle for each residue
        RD: numpy array (2x|R|) containing the surface depth evaluated from
            MSMS
        Surf: Surface from MSMS
        psaiaf: Features calculated from PSAIA (dictionary)
            'casa': chain asa
            'rasa': residue asa (total, backbone, sidechain, polar, non-polar)
            'rdpx': residue depth index
            'rhph': residue hydrophobicity index
            'rrasa': residue rasa
            'rcx': protrusion index
        resi: dictionary object (chainID,residueId)-> index in R.
            Here residueId (string) is residue number concatenated with
            its insertion code. This is used in the visaulization program.
    Notes:
        DSSP is used to calculate (SS,ASA,rASA,Phi,Psi).
        'nan' occurs whenever the value of a property is unavailable from DSSP.
    """
    def __init__(self,fname,auxdir=None):
        """
        usage: L=myPDB(fname,name)
            fname: path to the file
            name: name of the structure (optional, default: file name)
            auxdir: auxilliary information directory (should contain the pssm files (.mat))
        Instantiates the object.
        """
        """
        Set parameters:
        """
        self.HWS_seqFV=10.0
        self.HWS_PSSM=5.0

        (dpath,name,ftype)=getFileParts(fname)
        if auxdir is None:
            auxdir=dpath
        self.ftype=ftype.lower()
        self.name=name
        self.ifname=fname
        self.resi=dict()
        if self.ftype == '.pdb':
            (self.stx,self.R,self.pp,self.seq,self.S2Ri)=readPDB(fname)
            self.FV=getSeqFV(self.stx,self.R,self.HWS_seqFV)
            self.Coords=getCoords(self.R)
            self.__getStride()
            self.getSimMtx()
            self.__getBvalues()
            self.__calcRD()
            self.__getBlosumFeats()
            for (idx,r) in enumerate(self.R):
                self.resi[getResiId(r.get_full_id())]=idx
            self.__assignAux(auxdir)
        elif self.ftype == '.fasta':
            self.seq=readFASTA(fname)
            self.R=self.seq
            self.S2Ri=range(len(self.R))
            self.__getPSSM(auxdir)
            for (idx,r) in enumerate(self.R):
                self.resi[('A',str(idx))]=idx
        else:
            estring = 'Error: Unknown File Type. Input file must be PDB or FASTA.'
            print estring
            raise Exception(estring)
    def __getBlosumFeats(self):
        """
        Constructs the blosum 62 matrix representation (windowed)
        """
        wbsmf=getSubMatFeats(self.seq,self.HWS_PSSM)
        self.wbsmf=np.zeros((wbsmf.shape[0],len(self.R)))
        for k in range(len(self.seq)):
            idx=self.S2Ri[k]
            self.wbsmf[:,idx]=wbsmf[:,k]

    def __assignAux(self,auxpath):
        """
        Assign any features that require external precomputed files
        Specifically we use the profile (.mat) file computed using blastpgp (processBLASTPGP.py)
        to produce ASA,RASA and PSSM features
        """
        try:
            self.__getPSAIA(auxpath)
        except Exception as e:
            print "Error getting psaia features for",self.name,":",e," Is the .tbl file present?"

            # TODO(RAPH) removed this
        try:
            self.__getPSSM(auxpath)
        except Exception as e:
            print "Error getting profile features for",self.name,":",e," Is the .mat file present?"
    def save2FASTA(self,fpath,saveHeader=False):
        f = open(fpath, 'w+')
        if saveHeader:
            f.write('>'+self.name+'\n')
        f.write(self.seq)
        f.close()
    def __getPSSM(self,auxpath):
        #       Assigning PSSM Profile and using spinex

        ppath=os.path.join(auxpath,self.name+'.mat')
        if not (os.path.exists(ppath) and os.path.isfile(ppath)): #If mat file existsimpor
            print 'Extracting PSI BLAST profile'
            if self.ftype == '.pdb':
                fpath=os.path.join(auxpath,self.name+'.fasta')
                self.save2FASTA(fpath)
            elif self.ftype == '.fasta':
                fpath=self.ifname
            runPSIBLAST(fpath,db='nr',ofile=auxpath+self.name,niter=3)
        else:
            print 'Using existing PSI-BLAST Profile',ppath
        #pdb.set_trace()
        pfile=parsePSSMfile(ppath)#parsed file
        # TODO(Raph) might need to renable this.
#        zspinex=spinex.pssm_to_spinex(ppath)
        N=len(self.R)

        #pdb.set_trace()
        if pfile is not None: # and zspinex is not None:
            (pssm,psfm,info)=pfile
            wpssm=getWPSSM(pssm,self.HWS_PSSM)
            wpsfm=getWPSSM(psfm,self.HWS_PSSM)
            self.pssm=np.zeros((20,N))
            self.psfm=np.zeros((20,N))
            self.wpssm=np.zeros((wpssm.shape[0],N))
            self.wpsfm=np.zeros((wpsfm.shape[0],N))
            self.asa=np.zeros(N);self.asa.fill(np.nan)
            self.rasa=np.zeros(N);self.rasa.fill(np.nan)
            self.info=np.zeros(N);self.info.fill(np.nan)
#            (asa,rasa,_,_,_)=zspinex
            for k in range(len(self.seq)):
                idx=self.S2Ri[k]
                self.pssm[:,idx]=pssm[:,k]
                self.psfm[:,idx]=psfm[:,k]
                self.wpssm[:,idx]=wpssm[:,k]
                self.wpsfm[:,idx]=wpsfm[:,k]
                self.info[idx]=info[k]
#                self.asa[idx]=asa[k]
#                self.rasa[idx]=rasa[k]
        else:
            print 'Encountered some error processing spinex or psi-blast',self.name
            raise Exception(('Encountered some error processing spinex or psi-blast',self.name))

    def getSpectralVector(self):
        C=[]
        for r in self.R:
            C.append(r['CA'].get_coord())
        C=np.array(C)
        D=squareform(pdist(C))
        S=np.exp(-0.05*D)
        e_val,e_vec = LA.eig(S)
        v=(e_vec[:,0])
        dv=np.abs(v-v[:,np.newaxis])
        return (S,dv)
    def save(self,ofname=None,bdir='./'):
        """
        Save the object
        """
        if ofname is None:
            ofname=bdir+self.name+'.pdb.dill'
        output = open(ofname, 'wb')
        dill.dump(self, output,-1)
        output.close()

    @classmethod
    def loader(self,pklfile):
        """
        Load the class object from a pickel file
        """
        return cPickle.load(open(pklfile, "rb" ) )

    def getSimMtx(self,sg=2.0,thr=1e-3):
        """
        Instantiates the distance based similarity matrix (S). S is a tuple of
        lists (I,V). |I|=|V|=|R|. Each I[r] refers to the indices
        of residues in R which are "close" residue indexed by r in R and V[r]
        contains a list of the similarity scores for the corresponding residues.
        The distance between two residues is defined to be the minimum distance of
        any of their atoms. The similarity score is evaluated as
            s = exp(-d^2/(2*sg^2))
        This ensures that the range of similarity values is 0-1. sg (sigma)
        determines the extent of the neighborhood.
        Two residues are defined to be close to one another if their similarity
        score is greater than a threshold (thr).
        Residues (or ligands) for which DSSP features are not available are not
        included in the distance calculations.
        """
        self.sg=sg
        self.thr=thr
        sg=2*(sg**2)
        I=[[] for k in range(len(self.Coords))]
        V=[[] for k in range(len(self.Coords))]
        for i in range(len(self.Coords)):
            for j in range(i,len(self.Coords)):
                d=spatial.distance.cdist(self.Coords[i], self.Coords[j]).min()
                s=np.exp(-(d**2)/sg)
                if(s>thr):# and not np.isnan(self.Phi[i]) and not np.isnan(self.Phi[j])
                    I[i].append(j)
                    V[i].append(s)
                    if i!=j:
                        I[j].append(i)
                        V[j].append(s)
        self.S=(I,V)
        self.CN=np.array([len(a) for a in self.S[0]])
        self.__getHSE()

    def __getHSE(self):
        """
        Compute the Half sphere exposure statistics
        The up direction is defined as the direction of the side chain and is
        calculated by taking average of the unit vectors to different side chain
        atoms from the C-alpha atom
        Anything within the up half sphere is counted as up and the rest as
        down
        """
        N=len(self.R)
        Na=len(AA)
        self.UN=np.zeros(N)
        self.DN=np.zeros(N)
        self.UC=np.zeros((Na,N))
        self.DC=np.zeros((Na,N))
        for (i,r) in enumerate(self.R):
            u=getSideChainV(r)
            if u is None:
                self.UN[i]=np.nan
                self.DN[i]=np.nan
                self.UC[:,i]=np.nan
                self.DC[:,i]=np.nan
            else:
                idx=aaidx[getResLetter(r)]
                self.UC[idx,i]=self.UC[idx,i]+1
                self.DC[idx,i]=self.DC[idx,i]+1
                n=self.S[0][i]
                for j in n:
                    r2=self.R[j]
                    if is_aa(r2) and r2.has_id('CA'):
                        v2=r2['CA'].get_vector()
                        scode=getResLetter(r2)
                        idx=aaidx[scode]
                        if u[1].angle((v2-u[0])) < np.pi/2.0:
                            self.UN[i]=self.UN[i]+1
                            self.UC[idx,i]=self.UC[idx,i]+1
                        else:
                            self.DN[i]=self.DN[i]+1
                            self.DC[idx,i]=self.DC[idx,i]+1
        self.UC=self.UC/(1.0+self.UN)
        self.DC=self.DC/(1.0+self.DN)

    def __getDSSP(self):
        """
        Use DSSP to compute SS,ASA,rASA,Phi,Psi
        """
        ssk="HBEGITS-" #list of possible secondary structures
        sskey=dict(zip(ssk,range(len(ssk))))
        dssp=DSSP(self.stx[0],self.ifname)
        N=len(self.R)
        self.SS=np.zeros((len(ssk),N))
        self.ASA=np.zeros(N)
        self.rASA=np.zeros(N)
        self.Phi=np.zeros(N)
        self.Psi=np.zeros(N)
        for (idx,r) in enumerate(self.R):
            (_,_,cid,rid)=r.get_full_id()
            key=(cid,rid)
            if dssp.has_key(key):
                (_,s,self.ASA[idx],self.rASA[idx],self.Phi[idx],self.Psi[idx])=dssp[key]
                self.SS[sskey[s],idx]=1
            else:
                #print('here')
                #pdb.set_trace()
                self.SS[:,idx]=np.nan
                self.ASA[idx]=np.nan
                self.rASA[idx]=np.nan
                self.Phi[idx]=np.nan
                self.Psi[idx]=np.nan

    def __getStride(self):
        #ssk="HBEGITS-" #list of possible secondary structures
        ssk="HGIEBTC"
        sskey=dict(zip(ssk,range(len(ssk))))
        stridex=stride.stride_dict_from_pdb_file(self.ifname)
        N=len(self.R)
        self.SS=np.zeros((len(ssk),N))
        self.ASA=np.zeros(N)
        self.rASA=np.zeros(N)
        self.Phi=np.zeros(N)
        self.Psi=np.zeros(N)
        for (idx,r) in enumerate(self.R):
            key=getResiId(r.get_full_id())
            #pdb.set_trace()
            if stridex.has_key(key):
                # (aa,ss,phi,psi,asa,rasa)
                (_,s,self.Phi[idx],self.Psi[idx],self.ASA[idx],self.rASA[idx])=stridex[key]
                if not sskey.has_key(s):
                    print "unknown scondary structure! Add to dictionary!"
                    pdb.set_trace()
                self.SS[sskey[s],idx]=1
            else:
                print('key not found in stride processing!')
                #pdb.set_trace()
                self.SS[:,idx]=np.nan
                self.ASA[idx]=np.nan
                self.rASA[idx]=np.nan
                self.Phi[idx]=np.nan
                self.Psi[idx]=np.nan
    def __getPSAIA(self,auxpath=''):
        #fname=glob.glob(os.path.join(auxpath,self.name+'*bound.tbl'))[0]
        #pdict=psaia.make_psaia_dict(fname)
        pdict=psaia.runPSAIA(self.ifname)
        N=len(self.R)
        self.psaiaf={}
        self.psaiaf['casa']=np.zeros((5,N));self.psaiaf['casa'].fill(np.nan);
        self.psaiaf['rasa']=np.zeros((5,N));self.psaiaf['rasa'].fill(np.nan);
        self.psaiaf['rrasa']=np.zeros((5,N));self.psaiaf['rrasa'].fill(np.nan);
        self.psaiaf['rdpx']=np.zeros((6,N));self.psaiaf['rdpx'].fill(np.nan);
        self.psaiaf['rcx']=np.zeros((6,N));self.psaiaf['rcx'].fill(np.nan);
        self.psaiaf['rhph']=np.zeros(N);self.psaiaf['rhph'].fill(np.nan);
        for (idx,r) in enumerate(self.R):
            key=getResiId(r.get_full_id())
            #pdb.set_trace()
            if pdict.has_key(key):
                (casa,rasa,rrasa,rdpx,rcx,rhph)=pdict[key]
                self.psaiaf['casa'][:,idx]=casa
                self.psaiaf['rasa'][:,idx]=rasa
                self.psaiaf['rrasa'][:,idx]=rrasa
                self.psaiaf['rdpx'][:,idx]=rdpx
                self.psaiaf['rcx'][:,idx]=rcx
                self.psaiaf['rhph'][idx]=rhph
            else:
                print('key not found in PSAIA processing!')


    def __getBvalues(self):
        """
        Save the b value for each residue
        """
        self.B=np.zeros(len(self.R))
        for (idx,r) in enumerate(self.R):
            #pdb.set_trace()
            self.B[idx]=max([ak.get_bfactor() for ak in r])

    def __calcRD(self):
        """
        Get residue depth and surface
        """

        try:
            rd=ResidueDepth(self.stx[0])
        except:
            e = sys.exc_info()[0]
            print e
            #pdb.set_trace()
            print  "Oops!  MSMS not found.  Install MSMS and try again..."
            print  "Currently not using residue depth and surface features."
            return
        self.RD=np.empty((2,len(self.R)))
        self.RD.fill(np.nan)
        for (idx,r) in enumerate(self.R):
            (_,_,c,(h,rn,ic))=r.get_full_id()
            kk=(c,(h,rn,ic))
            if rd.has_key(kk):
                rdv=rd[kk]
                self.RD[0,idx]=rdv[0]
                self.RD[1,idx]=rdv[1]
        #self.rd=np.array(self.rd)
        """
        Get surface: It creates a single surface for all chains in the file
        """
        self.Surf=get_surface(self.stx[0])

    def showSurface(self):
        import pylab
        from mpl_toolkits.mplot3d import Axes3D
        fig = pylab.figure()
        ax = Axes3D(fig)
        ax.scatter(self.Surf[:,0],self.Surf[:,1],self.Surf[:,2],'.')
        plt.show()

    def plotRamachandran(self):
        """
        Plot the ramachandran plot (Psi vs. Phi)
        """
        plt.plot(self.Phi,self.Psi,'.')
        plt.xlabel("Phi")
        plt.ylabel("Psi")
        plt.show()


    def toFile(self,ofname=None):
        """
        saves the structure stx to the pdb file path (ofname)
        """
        if ofname is None:
            ofname=self.ifname
        io=PDBIO()
        io.set_structure(self.stx)
        io.save(ofname)

    def __str__( self ):
        s="myPDB object with Original File Name: "+self.name
        return s

if __name__=="__main__":
    usage="python myPDB.py ifile profile_path [ofile.pdb.pkl]"
    ofile=None
    alist=sys.argv
    if len(alist)==4:
        ofile=alist[3]
    if len(alist)<3:
        print 'Error: Improper arguments. Proper Usage is: ', usage
        exit(1)
    ifile=alist[1]
    #pdbpath='./pdb/'#'/media/sf_Desktop/Cam_test/'#'/s/chopin/b/grad/minhas/Desktop/PIANOX/benchmark4/structures/Shandar/'#
    #fid='1ATN_l_u.fasta'
    ppath=alist[2]
    L=myPDB(ifile,ppath)
    L.save(ofile)

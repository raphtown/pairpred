# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 21:55:38 2012
@author: Afsar
This module constructs positive and negative training examples from a data set
of pdb files 
"""
import glob
import os
import cPickle
from getPair import *
import itertools
import random
import numpy as np
class getExamplesDBD:
    """
    An object of this class will contain two dictionaries Pex amd Nex with 
    the same keys (complex_ids)
    Pex[key][0] is a list of tuples (L.R_idx,R.R_idx) of indicies from L and R
        proteins for complex id key that lie within a distance of 'dthr' of one 
        another. 
    Pex[key][1] is the number of residues in L
    Pex[key][2] is the number of residues in R
    Nex[key] is a list of tuples of indices from L and R that constitute negative
        examples
    dir: is the directory in which the pdb files were located from which the
        object has been constructed
    dthr: distance threshold
    Note: If there is an error in processing any file in a complex, the complex
        is not included in the object
    """
    def __init__(self,idir,dthr=6.0):
        """
        idir: path to the directory in which the files are located
        dthr: distance threshold (default: 6.0 Angstroms)
        """
        self.dthr=dthr
        self.dir=idir
        ufiles=glob.glob(os.path.join(idir,"*_l_u.pdb"))
        self.Pex=dict()
        self.Nex=dict()
        self.Plr=dict()
        for f in ufiles:            
            (_,cname,_)=getFileParts(f)
            cname=cname[:-4]
            print "############################ Processing: "+cname
            try:
                luname=os.path.join(idir,cname+'_l_u.pdb')
                runame=os.path.join(idir,cname+'_r_u.pdb')
                lbname=os.path.join(idir,cname+'_l_b.pdb')
                rbname=os.path.join(idir,cname+'_r_b.pdb')
                C=myPDBComplex([lbname,rbname],dthr)
                (D,LenU)=C.findNforU([luname,runame])            
                D=D[0][1]
                lenL=LenU[0]
                lenR=LenU[1]
                P=[]
                Pl=set([])
                Pr=set([])
                for (i,j,d) in D:
                    if ~ (np.isnan(i) or np.isnan(j)):
                        Pl.add(i)
                        Pr.add(j)
                        P.append((i,j))                
            except:
                print "!!!!!!!!!!!!!!!!!!!!!!!! Error processing: "+cname
                continue
            self.Nex[cname]=[]
            self.Pex[cname]=(P,lenL,lenR)
            self.Plr[cname]=[Pl,Pr]
    def save(self,ofname=None):
        """
        Save the object
        """
        if ofname is None:
            ofname='E_'+str(self.dthr)+'.lbl.pkl'
        print "File Saved: "+ofname
        output = open(ofname, 'wb')
        cPickle.dump(self, output)
        output.close()
        
    @classmethod   
    def loader(self,pklfile):
        """
        Load the class object from a pickel file
        """
        return cPickle.load(open(pklfile, "rb" ) )
        
    def getNegEx(self,cname):
        """
        Returns all negative examples for a given complex 'cname'
        LxR - P
        """
        try:
            L=[range(self.Pex[cname][1]),range(self.Pex[cname][2])]
            X=set(itertools.product(*L))
            return X.difference(set(self.Pex[cname][0]))
        except:
            return None
            
    def selectNegEx(self,cname=None,p=100.0,M=1000):    
        """
        Samples negative examples
        cname: list of complex ids for which the sampling is to be done
            if not specified, then the function automatically performs the
            sampling for all complexes in the object
        M: selection of sampling technique
            if M is not specified, then M=1000 and the minimum of M and 
                the number of negative examples in a complex is taken
            if M is set to None, then all negative examples are included
            if M is '+' then the number of negative examples is taken to 
                be equal to the minimum of the p (default 100) % of the number of positive examples and 
                the number of all possible negative examples
            Otherwise, M is expected to be a number and the number of negative
                examples is taken to be the minimum of p% of the total number
                of negative examples and M
            Negative examples are then randomly sampled and stored in Nex
        """
        if cname is None:
            cname=self.Pex.keys()
        if M is None:
            M=np.inf            
        for c in cname:            
            if M=='+-':
                Pl,Pr=self.Plr[c] #Positive examples in L and R
                Nl=range(self.Pex[c][1]) #All examples in L
                Nr=range(self.Pex[c][2]) #All examples in R
                XPlPr=set(itertools.product(Pl,Pr)).difference(set(self.Pex[c][0])) #Negative examples from (+,+) pairs
                XPlNr=(set(itertools.product(Pl,Nr)).difference(set(self.Pex[c][0]))).difference(XPlPr) #Negative examples from (+,-) pairs
                XNlPr=(set(itertools.product(Nl,Pr)).difference(set(self.Pex[c][0]))).difference(XPlPr) #Negative examples from (-,+) pairs
                XPN=XPlNr.union(XNlPr) #negative examples involving exactly 1 positive and 1 negative
                XNlNr=set(itertools.product(Nl,Pr)).difference(set(self.Pex[c][0])) #Negative examples from (-,-)
                XNlNr=XNlNr.difference(XPlPr)
                XNlNr=XNlNr.difference(XPN)
                npx=int(p*len(self.Pex[c][0])/100.0)
                #(+,+): 30%, (+,-):40%, (-,-): 30%
                npp=int(min(0.2*npx,len(XPlPr)))
                N=set(random.sample(XPlPr,npp))
                npn=int(min(0.3*npx,len(XPN)))
                N=N.union(set(random.sample(XPN,npn)))
                nn=npx-len(N)
                N=N.union(set(random.sample(XPN,nn)))
                self.Nex[c]=list(N)                        
            else:
                N=self.getNegEx(c)
                if M=='+':
                    n=int(min(p*len(self.Pex[c][0])/100.0,len(N)))
                else:
                    n=int(min(np.floor(p*len(N)/100.0),M))
                self.Nex[c]=random.sample(N,n)
if __name__=="__main__":    
    pdbpath='/s/chopin/b/grad/minhas/Desktop/DBD4N/DBD4/'
    E=getExamplesDBD(pdbpath)
    E.selectNegEx(p=150.0,M='+-') 
    E.save('EPN6_200.lbl.pkl')
            
            
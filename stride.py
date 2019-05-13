# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:50:27 2012

@author: Afsar
"""
import numpy as np
import tempfile
import os
from Bio.PDB import *
def getMAXASA(s=None):
    """
    This function returns a dictionary containing the maximum ASA for 
    different residues. when s=single, single letter codes of aa are also
    added to the dictionary
    """
    MAX_ACC={}
    MAX_ACC["ALA"]=106.0
    MAX_ACC["CYS"]=135.0
    MAX_ACC["ASP"]=163.0
    MAX_ACC["GLU"]=194.0
    MAX_ACC["PHE"]=197.0
    MAX_ACC["GLY"]=84.0
    MAX_ACC["HIS"]=184.0
    MAX_ACC["ILE"]=169.0
    MAX_ACC["LYS"]=205.0
    MAX_ACC["LEU"]=164.0
    MAX_ACC["MET"]=188.0
    MAX_ACC["ASN"]=157.0
    MAX_ACC["PRO"]=136.0
    MAX_ACC["GLN"]=198.0
    MAX_ACC["ARG"]=248.0
    MAX_ACC["SER"]=130.0
    MAX_ACC["THR"]=142.0
    MAX_ACC["VAL"]=142.0
    MAX_ACC["TRP"]=227.0
    MAX_ACC["TYR"]=222.0
    if s is not None and s is 'single':                
        for k in MAX_ACC.keys():
            MAX_ACC[to_one_letter_code[k]]=MAX_ACC[k]
    return MAX_ACC
def stride_dict_from_pdb_file(in_file, STRIDE="stride"):
    """
    Create a Stride dictionary from a PDB file.

    Example:
        >>> stride_dict=stride_dict_from_pdb_file("1fat.pdb")
        >>> (aa,ss,phi,psi,asa,rasa)=stride_dict[('A', 1)]

    @param in_file: pdb file
    @type in_file: string

    @param STRIDE: stride executable (argument to os.system)
    @type STRIDE: string

    @return: a dictionary that maps (chainid, resid) to 
        (aa,ss,phi,psi,asa,rasa)
    @rtype: {}
    #EXample: 
        {('A', '1'): ('GLY', 'C', 360.0, 119.38, 128.2, 1.0),
         ('A', '10'): ('ILE', 'E', -115.8, 136.5, 0.0, 0.0),...}
    Secondary structure codes:
        H	    Alpha helix
	  G	    3-10 helix
	  I	    PI-helix
	  E	    Extended conformation
	  B or	b   Isolated bridge
	  T	    Turn
	  C	    Coil (none of the above)
   IMPORTANT NOTE: if the protein chain	identifier is '	' (space), it
	   will	be substituted by '-' (dash) everywhere	in the STRIDE output.
	   The same is true  for  command  line	 parameters  involving	chain
	   identifiers where you have to specify '-' instead of	' '.
    """
    #import os

        
    def make_stride_dict(filename):

        """
        Return a stride dictionary that maps (chainid, resname, resid) to
        aa, ss and accessibility, from a stride output file.
        @param filename: the stride output file
        @type filename: string
        """
        MAX_ACC=getMAXASA()
        stride = {}
        handle = open(filename, "r")
        #print "herxxxxxxxxxxxxxxxxxxxxxxxxxxe"
        try:
            #kk=0
            for l in handle.readlines():
                #kk=kk+1
                #print kk
                sl = l.split()
                if sl[0] != "ASG": #if not detailed secondary structure record
                    continue
                #REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|  |-Area-|      ~~~~
                #ASG  ALA A    1    1    C          Coil    360.00    -35.26     120.7      ~~~~
                #0      1 2    3    4    5           6       7          8         9          10        
                # In cases where stride cannot recognize the residue type, it puts a '-' there
                # However, Bio.PDB uses ' ' so convert between the two                
                if sl[2]=='-':
                    sl[2]=' '
                
                resid=(sl[2],sl[3])
                aa=sl[1]
                ss=sl[5].upper() #There was b and B both from Bridge
                phi=float(sl[7])
                psi=float(sl[8])
                asa=float(sl[9])
                try:
                    rasa=asa/MAX_ACC[aa]
                    if rasa > 1.0: # we do get values greater than 1
                        rasa=1.0
                except KeyError:
                    rasa=np.nan
                stride[resid]=(aa,ss,phi,psi,asa,rasa)
                #construct a  key,value pair
                #pdb.set_trace()                
        finally:
            handle.close()
        return stride
        #return dssp, keys
    out_file = tempfile.NamedTemporaryFile(suffix='.stride')
    #out_file.flush()    # needed?
    out_file.close()
    ###########ADDED BY ME TO GET RID OF ACCESS DENIED ERROR   
        
    tname=out_file.name
    os.system("%s %s > %s" % (STRIDE, in_file, tname))
    out_dict= make_stride_dict(tname)
    os.remove(tname)
    return out_dict
    

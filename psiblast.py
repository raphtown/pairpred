# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 14:00:03 2013

@author: root
Description: Code for producing the profiles from a FASTA file using PSIBLAST
Dependencies: Installed psiblast (version 2.2.27)
    .ncbirc file in the home directory of the user containing the path to the databases
    Example of .ncbirc file
        [BLAST]
        BLASTDB=/s/chopin/c/proj/protfun/arch/x86_64/blast+/bin/nr
    The nr database was downloaded from ftp://ftp.ncbi.nih.gov/blast/db/. It was
    then untared and put in the above folder.
"""
import os
import pdb
import numpy as np
BLAST_DIR='/s/chopin/c/proj/protfun/arch/x86_64/ncbi-blast-2.2.27+-x64-linux/ncbi-blast-2.2.27+/bin/'#'/s/chopin/c/proj/protfun/arch/x86_64/blast+/bin/'
if BLAST_DIR not in  os.environ["PATH"]:
    print 'Adding '+BLAST_DIR+' to system path'
    os.environ["PATH"]+=os.pathsep+BLAST_DIR

def runPSIBLAST(f,db='nr',ofile=None,niter=3):
    blastdbpath='/scratch/PI/rondror/ppi/non-redundant/nr'
    if ofile is None:
        ofile = f
    #cmdstr='blastpgp -d nr -i '+f+' -Q '+f+'.mat -a 4 -j 3' #for version 2.2.24 (Aug-08-2010)
    #TODO(Raph): Log that you added blastdbpath here.
    cmdstr='psiblast -num_threads 4 -query '+f+' -db '+blastdbpath+' -out '+ofile+'.psi.txt'+' -num_iterations '+str(niter)+' -out_pssm '+ofile+'.pssm'+' -out_ascii_pssm '+ofile+'.mat_' #version 2.2.24+
    print cmdstr
    if os.system(cmdstr) == 0:
        print('PSIBLAST Successful : ', cmdstr)
        correctFile(ofile+'.mat_',ofile+'.mat')
        os.remove(ofile+'.mat_')
        return 0
    else:
        print('PSIBLAST Processing Failed.')
        return -1
def correctFile(ifile,ofile):
    """
    Blast 2.2.27 ignores the spacing in some lines
    """
    of=open(ofile,'w')
    for f in open(ifile,'r'):
        fs=f.split()
        if len(fs) and fs[0].isdigit():
            if len(fs)<44:
                fpssm=' '.join(split_str_into_len(f[9:69],l=3))
                f=f[:9]+fpssm+f[69:]
        of.write(f)
    of.close()
                #print f
                #pdb.set_trace()
def split_str_into_len(s, l=3):
    """ Split a string into chunks of length l """
    return [(s[i:i+l]).strip() for i in range(0, len(s), l)]

def parsePSSMfile(fname):
        pssm=[]
        psfm=[]
        info=[]
        try:
            for f in open(fname,'r'):
                f=f.split()
                if len(f) and f[0].isdigit(): #the first character must be a position
                    _pssm=[float(i) for i in f[2:22]]#split_str_into_len(f[9:69],l=3)#
                    pssm.append(_pssm)
                    #f=f[69:].split()
                    _psfm=[float(i)/100.0 for i in f[22:42]]#[float(i)/100.0 for i in f[:20]]#
                    psfm.append(_psfm)
                    #info.append(float(f[20]))
                    info.append(float(f[42]))
            z=(np.array(pssm).T,np.array(psfm).T,np.array(info))
        except IOError as e:
            print e
            z=None
        return z

if __name__=="__main__":
    f='1URSA'
    if not runPSIBLAST(f,db='nr'):
        z=parsePSSMfile(f+'.mat')
    #correctFile('2HRK_r_u.mat','temp.mat')

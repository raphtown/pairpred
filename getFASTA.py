# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 09:11:25 2012

@author: root
Description: Saves the sequence FASTA files from the structure files
"""

from BISEPutils import readPDB
from getExamplesDBD import *
import matplotlib.pyplot as plt
import pdb
try:
    import cPickle as pickle
except:
    import pickle
    
"""
This script saves the sequences from all pdb.pkl files in a folder (using the getExamplesDBD file)
"""

def saveFASTA(odir,luname):   
    """
    Given a pdbpkl file it saves the sequence of that into a fasta file with the same
    name as the pdbpkl file and no extension in a given folder
    """
    ofname=os.path.join(odir,getFileParts(luname)[1])
    (_,_,_,seq,_)=readPDB(luname)
    f = open(ofname, 'w+')
    f.write(seq)
    f.close()
    
#--------------------MAIN--------------------------------    
if __name__=="__main__":
    pdbpklpath='./Benchmark4'
    exfname=pdbpklpath+'/EB4_6.0.lbl.pkl'
    odir='./FASTA'
    E=getExamplesDBD.loader(exfname)
#    rigid=['1AHW', '1BVK', '1DQJ', '1E6J', '1JPS', '1MLC', '1VFB', '1WEJ', '2FD6', '2I25', '2VIS', '1BJ1', '1FSK', '1I9R', '1IQD', '1K4C', '1KXQ', '1NCA', '1NSN', '1QFW', '1QFW', '2JEL', '1AVX', '1AY7', '1BVN', '1CGI', '1CLV', '1D6R', '1DFJ', '1E6E', '1EAW', '1EWY', '1EZU', '1F34', '1FLE', '1GL1', '1GXD', '1HIA', '1JTG', '1MAH', '1N8O', '1OC0', '1OPH', '1OYV', '1OYV', '1PPE', '1R0R', '1TMQ', '1UDI', '1YVB', '2ABZ', '2B42', '2J0T', '2MTA', '2O8V', '2OUL', '2PCC', '2SIC', '2SNI', '2UUY', '3SGQ', '4CPA', '7CEI', '1A2K', '1AK4', '1AKJ', '1AZS', '1B6C', '1BUH', '1E96', '1EFN', '1F51', '1FC2', '1FCC', '1FFW', '1FQJ', '1GCQ', '1GHQ', '1GLA', '1GPW', '1H9D', '1HCF', '1HE1', '1I4D', '1J2J', '1JWH', '1K74', '1KAC', '1KLU', '1KTZ', '1KXP', '1ML0', '1OFU', '1PVH', '1QA9', '1RLB', '1RV6', '1S1Q', '1SBB', '1T6B', '1US7', '1WDW', '1XD3', '1XU1', '1Z0K', '1Z5Y', '1ZHH', '1ZHI', '2A5T', '2A9K', '2AJF', '2AYO', '2B4J', '2BTF', '2FJU', '2G77', '2HLE', '2HQS', '2OOB', '2OOR', '2VDB', '3BP8', '3D5S']
#    medium=['1BGX', '1ACB', '1IJK', '1JIW', '1KKL', '1M10', '1NW9', '1GP2', '1GRN', '1HE8', '1I2M', '1IB1', '1K5D', '1LFD', '1MQ8', '1N2C', '1R6Q', '1SYX', '1WQ1', '1XQS', '1ZM4', '2CFH', '2H7V', '2HRK', '2J7P', '2NZ8', '2OZA', '2Z0E', '3CPH']
#    hard=['1E4K', '2HMI', '1F6M', '1FQ1', '1PXV', '1ZLI', '2O3B', '1ATN', '1BKD', '1DE4', '1EER', '1FAK', '1H1V', '1IBR', '1IRA', '1JK9', '1JMO', '1JZD', '1R8S', '1Y64', '2C0L', '2I9B', '2IDO', '2OT3']
    cids=E.Pex.keys()
    for cid in cids:
        print "Processing: ",cid
        try:        
            saveFASTA(odir,os.path.join(pdbpklpath,cid+'_l_u.pdb'))
            saveFASTA(odir,os.path.join(pdbpklpath,cid+'_r_u.pdb'))
        except Exception as e:
            print e
            continue
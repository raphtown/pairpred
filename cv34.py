# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 21:27:13 2013
Description: To evaluate results when training is  performed on DBD 3 and testing on DBD4 and vice versa
@author: root
"""

from dbdscrpp3 import *
if __name__=="__main__":
    pwKfname='dbdKernel.dbk.pkl'    
    f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN', '1N2C', '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
    f4=['2A5T', '3CPH', '1ZHH', '2ABZ', '1LFD', '2OUL', '1JIW', '2B4J', '1SYX', '1FLE', '1JTG', '2AYO', '4CPA', '1CLV', '1OC0', '1XU1', '1R6Q', '2O3B', '1US7', '3D5S', '1JZD', '1HCF', '1OYV', '2OZA', '1H9D', '2A9K', '2J0T', '2Z0E', '3BP8', '2IDO', '1WDW', '1ZLI', '2VDB', '1RV6', '1FFW', '1F6M', 'BOYV', '1JWH', '2OOR', '1MQ8', '1GL1', '1PVH', '2I9B', '1OFU', '1GXD', '3SGQ', '1JK9', '1ZM4', '1FCC', '2G77', '2J7P', '2FJU']
    #f3=['1XQS','1A2K']
    #f4=['1ACB','1AHW','1XXX'] #'1A2K','1ACB','1AHW','1XQS
    print "Loading Kernel from file",pwKfname
    s=SVM()
    s.C=10  
    t1=time()
    mK=getKernelData(*dbdKernel.loader(pwKfname)) #NOTE THAT OBJECT IS BEING LOADED FROM FILE
    t2=time()
    print "Pairwise kernel loaded in",t2-t1,"s"   
    #pdb.set_trace()
    ttP,dkey=getComplexFoldList(mK,nfolds=[f3,f4])
    mK.attachLabels(Labels(mK.labels.L))    #pyml did not work with the string labels
    print "Starting CV ... "
    r=mycvFromFolds(s,mK,testingPatterns=ttP)
    if r is not None:
        output = open('result.res.pkl', 'wb')
        cPickle.dump((r,dkey), output)
        output.close()
        (auc,(fp,tp),(A,Rx,D,L,cids,r,dkey))=getAUC((r,dkey))        
        print "CV Complete. Total Time taken (s):",time()-t1
        print "AUC for testing on DBD 3:",np.mean([A[i] for i,c in enumerate(cids) if c in f3])
        print "AUC for testing on DBD 4:",np.mean([A[i] for i,c in enumerate(cids) if c in f4])
        print "Complex-wise AUC = ",auc
        print "Overall AUC = ", r.getROC()
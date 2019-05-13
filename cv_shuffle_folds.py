# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 19:24:11 2013

@author: root
"""

from dbdscrpp3 import *

if __name__=="__main__":      
    pwKfname='869dbdKernel.dbk.pkl'
    print "Loading Kernel from file ",pwKfname
    t1=time()
    mK=getKernelData(*dbdKernel.loader(pwKfname)) #NOTE THAT OBJECT IS BEING LOADED FROM FILE
    t2=time()
    print "Pairwise kernel loaded in",t2-t1,"s"
    Nshuffles=5
    ttP=[[] for n in range(Nshuffles)]
    dkey=[[] for n in range(Nshuffles)]
    for n in range(Nshuffles):
        ttP[n],dkey[n]=getComplexFoldList(mK,shuffle=True)

    mK.attachLabels(Labels(mK.labels.L))    #pyml did not work with the string labels
    aucx=[[] for n in range(Nshuffles)]
    for n in range(Nshuffles):
        print "Starting CV for shuffle ",n
        s=SVM()
        s.C=10 
        r=mycvFromFolds(s,mK,testingPatterns=ttP[n])
        fname='result_'+str(n)+'.res.pkl'
        if r is not None:
            output = open(fname, 'wb')
            cPickle.dump((r,dkey[n]), output)
            output.close()
            (auc,_,_)=getAUC((r,dkey[n]))
            aucx[n]=auc
            print "CV Complete for shuffle",n,". Total Time taken (s):",time()-t1
            print "Complex-wise average AUC",auc
            print "Overall AUC = ", r.getROC()
    if r is not None:
        print 'MEAN AUC =',np.mean(aucx)
        print 'STD AUC =',np.std(aucx)
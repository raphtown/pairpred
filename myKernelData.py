# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 16:02:03 2012

@author: root
"""
from PyML import *
from PyML.containers.ext import ckernel
from PyML.containers.ext import ckerneldata
from PyML.containers.baseDatasets import WrapperDataSet
from PyML.containers import labels
import numpy as np
import pdb
import cPickle
class myKernelData(WrapperDataSet, ckerneldata.KernelData):
    """
    My implementation of
    /usr/local/lib/python2.7/dist-packages/PyML/containers/kernelData.py
    to handle direct conversion of a Kernel array from numpy into 
    KernelData style Object which can then be used for training / cv etc.
    Usage:
        mK=myKernelData(K,[labels=Y,patternID=S])
    """
    isVector = False
    def __init__(self, arg = None, **args) :
        """
        
        """        
        if arg.__class__ == self.__class__ :
            self.copyConstruct(arg, **args)			
        elif type(arg) == np.ndarray :
            self.constructFromMatrix(arg, **args)
            WrapperDataSet.attachKernel(self)
    def copy(self, other, patterns, deepcopy) :        
        ckerneldata.KernelData.__init__(self, other, patterns)	
        
    def constructFromMatrix(self, K, **args) :		
        matrix = ckernel.KernelMatrix()
        matrix.thisown = 0
        for idx in range(K.shape[0]):
            matrix.addRow(K[idx,:])
        ckerneldata.KernelData.__init__(self, matrix)        
        if 'labels' in args :
            self.attachLabels(labels.Labels(args['labels'], **args))
        else :
            patternID=range(K.shape[0])
            self.attachLabels(labels.Labels(None, patternID = patternID))
            
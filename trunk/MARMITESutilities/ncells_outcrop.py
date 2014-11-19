# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 22:16:20 2014

@author: apf
"""

import numpy as np

ar1 = np.asarray([[1,1,0], [0,1,0], [0,1,1], [1,0,0], [1,1,1]], dtype = np.int16)
ar2 = np.asarray([[1,1,0], [1,1,1], [0,1,0], [0,0,1], [1,1,1]], dtype = np.int16)
ar3 = np.asarray([[1,1,0], [1,1,1], [0,1,1], [0,1,1], [1,1,1]], dtype = np.int16)
ar4 = np.asarray([[1,0,1], [1,1,1], [1,1,0], [1,0,1], [1,1,1]], dtype = np.int16)

nlay = 4

ncell = []
ar = np.zeros(ar1.shape)
outcropL = np.zeros(ar1.shape)
for l, L in enumerate([ar1,ar2,ar3,ar4]):
    ar += L
    ncell.append((ar*L == 1).sum())
    outcropL += ((outcropL == 0) & (ar == 1))*(l+1)
    
print ncell
print outcropL
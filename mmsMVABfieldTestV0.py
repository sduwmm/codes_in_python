# -*- coding: utf-8 -*-
"""
Created on Tue May  21 18:26:28 2019

@author: Mengmeng Wang
"""

"""
Minimum Variance Analysis
Magnetic field
Verion Test
"""

"""
!!!Note

"""

import cdflib
import numpy as np
import math
import matplotlib.pyplot as plt
###
import time


pyStartTime = time.time()

timeInput1 = [2016,12,17,12,55, 4,000]
timeInput2 = [2016,12,17,12,55,43,000]

fileNameFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2016/12/17/mms1_fgm_brst_l2_20161217125504_v5.87.0.cdf')
infoFgm     = fileNameFgm.cdf_info()
epochFgm    = fileNameFgm.varget('Epoch')
bGSEVec     = fileNameFgm.varget('mms1_fgm_b_gse_brst_l2')

epochStartType2000 = cdflib.cdfepoch.compute_tt2000(timeInput1)
epochStopType2000  = cdflib.cdfepoch.compute_tt2000(timeInput2)
timeTuple1         = np.where(epochFgm > epochStartType2000)
timeTuple2         = np.where(epochFgm < epochStopType2000)
timeArr1           = timeTuple1[0]
timeArr2           = timeTuple2[0]
timeArrToCal       = np.intersect1d(timeArr1, timeArr2)

bxGSEVec = bGSEVec[:,0]
byGSEVec = bGSEVec[:,1]
bzGSEVec = bGSEVec[:,2]
btGSEVec = bGSEVec[:,3]

bxGSEToCal = bxGSEVec[timeArrToCal]
byGSEToCal = byGSEVec[timeArrToCal]
bzGSEToCal = bzGSEVec[timeArrToCal]
btGSEToCal = btGSEVec[timeArrToCal]

mm = np.ones((3, 3))
mm[0][0] = np.mean(bxGSEToCal * bxGSEToCal) - np.mean(bxGSEToCal) * np.mean(bxGSEToCal)
mm[0][1] = np.mean(bxGSEToCal * byGSEToCal) - np.mean(bxGSEToCal) * np.mean(byGSEToCal)
mm[0][2] = np.mean(bxGSEToCal * bzGSEToCal) - np.mean(bxGSEToCal) * np.mean(bzGSEToCal)
mm[1][0] = np.mean(byGSEToCal * bxGSEToCal) - np.mean(byGSEToCal) * np.mean(bxGSEToCal)
mm[1][1] = np.mean(byGSEToCal * byGSEToCal) - np.mean(byGSEToCal) * np.mean(byGSEToCal)
mm[1][2] = np.mean(byGSEToCal * bzGSEToCal) - np.mean(byGSEToCal) * np.mean(bzGSEToCal)
mm[2][0] = np.mean(bzGSEToCal * bxGSEToCal) - np.mean(bzGSEToCal) * np.mean(bxGSEToCal)
mm[2][1] = np.mean(bzGSEToCal * byGSEToCal) - np.mean(bzGSEToCal) * np.mean(byGSEToCal)
mm[2][2] = np.mean(bzGSEToCal * bzGSEToCal) - np.mean(bzGSEToCal) * np.mean(bzGSEToCal)
[D, V] = np.linalg.eig(mm)
#DD = np.ones((1, 3))
#DD[0][0] = D[0]
#DD[0][1] = D[1]
#DD[0][2] = D[2]
lambda1 = np.min(D)
Imin    = np.where(D == np.min(D, axis = 0))
lambda3 = np.max(D)
Imax    = np.where(D == np.max(D, axis = 0))
kLoop = np.arange(0,3)
for k in kLoop:
    if D[k] != lambda1 and D[k] != lambda3:
        lambda2 = D[k]
        Imid = k

Imin = Imin[0]
Imin = Imin[0]
Imax = Imax[0]
Imax = Imax[0]
N = V[:, Imin]
M = V[:, Imid]
L = V[:, Imax]

#def mva(bx, by, bz):
#    mm = np.ones((3,3))
#    mm[0][0] = np.mean(bxGSEToCal*bxGSEToCal) - np.mean(bxGSEToCal)*np.mean(bxGSEToCal)
#    mm[0][1] = np.mean(bxGSEToCal*byGSEToCal) - np.mean(bxGSEToCal)*np.mean(byGSEToCal)
#    mm[0][2] = np.mean(bxGSEToCal*bzGSEToCal) - np.mean(bxGSEToCal)*np.mean(bzGSEToCal)
#    mm[1][0] = np.mean(byGSEToCal*bxGSEToCal) - np.mean(byGSEToCal)*np.mean(bxGSEToCal)
#    mm[1][1] = np.mean(byGSEToCal*byGSEToCal) - np.mean(byGSEToCal)*np.mean(byGSEToCal)
#    mm[1][2] = np.mean(byGSEToCal*bzGSEToCal) - np.mean(byGSEToCal)*np.mean(bzGSEToCal)
#    mm[2][0] = np.mean(bzGSEToCal*bxGSEToCal) - np.mean(bzGSEToCal)*np.mean(bxGSEToCal)
#    mm[2][1] = np.mean(bzGSEToCal*byGSEToCal) - np.mean(bzGSEToCal)*np.mean(byGSEToCal)
#    mm[2][2] = np.mean(bzGSEToCal*bzGSEToCal) - np.mean(bzGSEToCal)*np.mean(bzGSEToCal)
#    [V,D] = np.linalg.eig(mm)
#    DD = np.ones((1,3))
#    DD[0][0] = D[0,0]
#    DD[1][1] = D[1,1]
#    DD[2][2] = D[2,2]
#    [lambda1, Imin] = min(DD)
#    [lambda3, Imax] = max(DD)
#    return(Imin,Imax)
#[v,d] = mva(bxGSEToCal, byGSEToCal, bzGSEToCal)
#print(v)
#print(d)

pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))

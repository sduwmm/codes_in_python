# -*- coding: utf-8 -*-
'''
Created on Tue May 27 19:45:30 2019

@author: Mengmeng Wang
'''

'''
Minimum Variance Analysis
Magnetic field
Version 0
'''

'''
!!!Note
'''

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
epochToCal = epochFgm[timeArrToCal]

def mva(bx,by,bz):
    mm       = np.ones((3, 3))
    mm[0][0] = np.mean(bx * bx) - np.mean(bx) * np.mean(bx)
    mm[0][1] = np.mean(bx * by) - np.mean(bx) * np.mean(by)
    mm[0][2] = np.mean(bx * bz) - np.mean(bx) * np.mean(bz)
    mm[1][0] = np.mean(by * bx) - np.mean(by) * np.mean(bx)
    mm[1][1] = np.mean(by * by) - np.mean(by) * np.mean(by)
    mm[1][2] = np.mean(by * bz) - np.mean(by) * np.mean(bz)
    mm[2][0] = np.mean(bz * bx) - np.mean(bz) * np.mean(bx)
    mm[2][1] = np.mean(bz * by) - np.mean(bz) * np.mean(by)
    mm[2][2] = np.mean(bz * bz) - np.mean(bz) * np.mean(bz)
    [D, V]  = np.linalg.eig(mm)
    lambda1 = np.min(D)
    Imin    = np.where(D == lambda1)
    lambda3 = np.max(D)
    Imax    = np.where(D == lambda3)
    kLoop   = np.arange(0,3)
    for k in kLoop:
        if D[k] != lambda1 and D[k] != lambda3:
            lambda2 = D[k]
            Imid    = k
    Imin = Imin[0]
    Imin = Imin[0]
    Imax = Imax[0]
    Imax = Imax[0]
    N    = V[:, Imin]
    M    = V[:, Imid]
    L    = V[:, Imax]
    return(lambda1, lambda2, lambda3, N, M, L, V)
[lambda1, lambda2, lambda3, N, M, L, arrRot] = mva(bxGSEToCal, byGSEToCal, bzGSEToCal)


bGSEToRot = np.column_stack((bxGSEToCal, byGSEToCal, bzGSEToCal))
bLMN      = np.dot(bGSEToRot, arrRot)

fig, ax = plt.subplots()
xlim    = [min(epochToCal), max(epochToCal)]
colors  = ['blue', 'green', 'red']
labels  = ['bn',   'bm',    'bl']
kLoop   = np.arange(0,3)
ax.set_xlim(xlim)
ax.set_ylim([-40,30])
for k in kLoop:
    ax.plot(epochToCal,bLMN[:,k], linewidth = 0.5, color = colors[k], label = labels[k])
    ax.legend(loc = 2, fontsize = 15, markerscale = 1.0, handlelength = 2, handleheight = 0.5)

xtick1 = [2016, 12, 17, 12, 55, 10, 000]
xtick2 = [2016, 12, 17, 12, 55, 20, 000]
xtick3 = [2016, 12, 17, 12, 55, 30, 000]
xtick4 = [2016, 12, 17, 12, 55, 40, 000]

xtick1 = cdflib.cdfepoch.compute_tt2000(xtick1)
xtick2 = cdflib.cdfepoch.compute_tt2000(xtick2)
xtick3 = cdflib.cdfepoch.compute_tt2000(xtick3)
xtick4 = cdflib.cdfepoch.compute_tt2000(xtick4)

epochTick = [xtick1, xtick2, xtick3, xtick4]
strTick = ['12:55:10', '12:55:20', '12:55:30', '12:55:40']
ax.set_xticks(epochTick)
ax.set_xticklabels(strTick, fontsize = 15)

pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))
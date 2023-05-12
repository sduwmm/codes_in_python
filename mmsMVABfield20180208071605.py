# -*- coding: utf-8 -*-
'''
Created on Tue May 27 19:45:30 2019

@author: Mengmeng Wang
'''

'''
Minimum Variance Analysis
Magnetic field
Version 0
case 20180208_071605
'''

'''
!!!Note
'''

import cdflib
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
###
import time

'''
function MVA:===========================================================================================================
'''
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

'''
!!!Note:
arrange eigenvalues from small to large
lambda1, lambda2, lambda3
corresponding to eigenvectors
N, M, L
corresponding to variation of magnetic field 
minimum, intermediate, maximum
!!!Note
'''
'''
function MVA:===========================================================================================================
'''

pyStartTime = time.time()

timeInput1 = [2018, 2, 8, 7, 16,  5,000]
timeInput2 = [2018, 2, 8, 7, 17, 15,000]

fileNameFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2018/02/08/mms1_fgm_brst_l2_20180208071603_v5.125.0.cdf')
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


[lambda1, lambda2, lambda3, N, M, L, arrRot] = mva(bxGSEToCal, byGSEToCal, bzGSEToCal)


bGSEToRot = np.column_stack((bxGSEToCal, byGSEToCal, bzGSEToCal))
bLMN      = np.dot(bGSEToRot, arrRot)

fig = plt.figure(figsize = (10, 8))
fig.subplots_adjust(top    = 0.94)
fig.subplots_adjust(bottom = 0.07)
fig.subplots_adjust(left   = 0.13)
fig.subplots_adjust(right  = 0.93)
gs  = GridSpec(3, 2, figure = fig)
ax1 = fig.add_subplot(gs[0,:])
ax2 = fig.add_subplot(gs[1,:])
ax3 = fig.add_subplot(gs[-1,0])
ax4 = fig.add_subplot(gs[-1,-1])
#plt.subplots_adjust(wspace = 0, hspace = 0)


xlim      = [min(epochToCal), max(epochToCal)]
colors    = ['blue', 'green', 'red']
labelsGSE = ['bx', 'by', 'bz']
labelsLMN = ['bn',   'bm',    'bl']
kLoop     = np.arange(0,3)

for k in kLoop:
    ax1.plot(epochToCal, bGSEToRot[:,k], linewidth = 0.8, color = colors[k], label = labelsGSE[k])
    ax1.legend(loc = 2, fontsize = 8, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
    ax2.plot(epochToCal,bLMN[:,k], linewidth = 0.8, color = colors[k], label = labelsLMN[k])
    ax2.legend(loc = 2, fontsize = 8, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
xtick1 = [2018, 2, 8, 7, 16, 15, 000]
xtick2 = [2018, 2, 8, 7, 16, 30, 000]
xtick3 = [2018, 2, 8, 7, 16, 45, 000]
xtick4 = [2018, 2, 8, 7, 17,  0, 000]

xtick1 = cdflib.cdfepoch.compute_tt2000(xtick1)
xtick2 = cdflib.cdfepoch.compute_tt2000(xtick2)
xtick3 = cdflib.cdfepoch.compute_tt2000(xtick3)
xtick4 = cdflib.cdfepoch.compute_tt2000(xtick4)

epochTick = [xtick1, xtick2, xtick3, xtick4]
strTick = ['04:16:15', '07:16:30', '07:16:45', '07:17:00']

ax1.set_xlim(xlim)
ax1.set_ylim([-35,35])
ax1.set_xticks(epochTick)
ax1.set_xticklabels(strTick, fontsize = 15)
ax2.set_xlim(xlim)
ax2.set_ylim([-40,30])
ax2.set_xticks(epochTick)
ax2.set_xticklabels(strTick, fontsize = 15)



timeSatrtHodo = [2018, 2, 8, 7, 16, 25, 000]
timeStopHodo  = [2018, 2, 8, 7, 16, 35, 000]
epochHodoStartType2000 = cdflib.cdfepoch.compute_tt2000(timeSatrtHodo)
epochHodoStopType2000  = cdflib.cdfepoch.compute_tt2000(timeStopHodo)
timeHodoTuple1         = np.where(epochToCal > epochHodoStartType2000)
timeHodoTuple2         = np.where(epochToCal < epochHodoStopType2000)
timeHodoArr1           = timeHodoTuple1[0]
timeHodoArr2           = timeHodoTuple2[0]
timeHodoArrToCal       = np.intersect1d(timeHodoArr1, timeHodoArr2)

bN = bLMN[:,0]
bM = bLMN[:,1]
bL = bLMN[:,2]
bNHodo = bN[timeHodoArrToCal]
bMHodo = bM[timeHodoArrToCal]
bLHodo = bL[timeHodoArrToCal]
ax3.plot(bMHodo, bLHodo, color = 'steelblue')
ax3.scatter(bMHodo[0], bLHodo[0], marker = 'o', color = 'steelblue')
ax3.scatter(bMHodo[bMHodo.size - 1], bLHodo[bLHodo.size - 1], marker = '>', color = 'steelblue')
ax3.set_xlabel('B(int)')
ax3.set_ylabel('B(max)')
ax4.plot(bNHodo, bLHodo, color = 'darkorange')
ax4.scatter(bNHodo[0], bLHodo[0], marker = 'o', linewidth = 5, color = 'darkorange')
ax4.scatter(bNHodo[bNHodo.size - 1], bLHodo[bLHodo.size -1], marker = '>', linewidth = 5, color = 'darkorange')
ax4.set_xlabel('B(min)')
ax4.set_ylabel('B(max)')



pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))
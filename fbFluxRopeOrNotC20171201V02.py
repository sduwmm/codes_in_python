# -*- coding: utf-8 -*-
'''
Created on Fri Feb 14 15:45:30 2020

@author: Mengmeng Wang
'''

'''
'''

'''
!!!Note
'''

import cdflib
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.gridspec import GridSpec
from numpy import linalg
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import time


'''
function MVA:===========================================================================================================
'''


def mva(bx, by, bz):
    mm = np.ones((3, 3))
    mm[0][0] = np.mean(bx * bx) - np.mean(bx) * np.mean(bx)
    mm[0][1] = np.mean(bx * by) - np.mean(bx) * np.mean(by)
    mm[0][2] = np.mean(bx * bz) - np.mean(bx) * np.mean(bz)
    mm[1][0] = np.mean(by * bx) - np.mean(by) * np.mean(bx)
    mm[1][1] = np.mean(by * by) - np.mean(by) * np.mean(by)
    mm[1][2] = np.mean(by * bz) - np.mean(by) * np.mean(bz)
    mm[2][0] = np.mean(bz * bx) - np.mean(bz) * np.mean(bx)
    mm[2][1] = np.mean(bz * by) - np.mean(bz) * np.mean(by)
    mm[2][2] = np.mean(bz * bz) - np.mean(bz) * np.mean(bz)
    [D, V] = np.linalg.eig(mm)
    lambda1 = np.min(D)
    Imin = np.where(D == lambda1)
    lambda3 = np.max(D)
    Imax = np.where(D == lambda3)
    kLoop = np.arange(0, 3)
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
    DNew = [D[Imin], D[Imid], D[Imax]]
    VNew = [N, M, L]
    return (lambda1, lambda2, lambda3, N, M, L, DNew, VNew)



pyStartTime = time.time()

figurePath = 'C:\\Users\\BRIAR\\Desktop\\OneDrive\\work\\code&figure\\'
if not os.path.exists(figurePath):
    os.makedirs(figurePath)

fig1 = plt.figure(figsize = (8, 10), dpi = 200)
fig1.subplots_adjust(top = 0.94)
fig1.subplots_adjust(bottom = 0.07)
fig1.subplots_adjust(left = 0.13)
fig1.subplots_adjust(right = 0.90)
gs = GridSpec(6,1, figure = fig1)
ax0 = fig1.add_subplot(gs[0, :])
ax1 = fig1.add_subplot(gs[1, :])
ax2 = fig1.add_subplot(gs[2, :])
ax3 = fig1.add_subplot(gs[3, :])
ax4 = fig1.add_subplot(gs[4, :])
ax5 = fig1.add_subplot(gs[5, :])


timeInput0 = [2017, 12, 1, 14, 40, 50, 000]
timeInput1 = [2017, 12, 1, 14, 41,  6, 000]



fileNameFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2017/12/01/mms1_fgm_brst_l2_20171201143933_v5.114.0.cdf')
infoFgm     = fileNameFgm.cdf_info()
epochFgm    = fileNameFgm.varget('Epoch')
bGSEVec     = fileNameFgm.varget('mms1_fgm_b_gse_brst_l2')

epochFgmStartType2000 = cdflib.cdfepoch.compute_tt2000(timeInput0)
epochFgmStopType2000  = cdflib.cdfepoch.compute_tt2000(timeInput1)
timeFgmTuple0         = np.where(epochFgm > epochFgmStartType2000)
timeFgmTuple1         = np.where(epochFgm < epochFgmStopType2000)
timeFgmArr0           = timeFgmTuple0[0]
timeFgmArr1           = timeFgmTuple1[0]
timeFgmArrToCal       = np.intersect1d(timeFgmArr0, timeFgmArr1)

bxGSEVec = bGSEVec[:,0]
byGSEVec = bGSEVec[:,1]
bzGSEVec = bGSEVec[:,2]
btGSEVec = bGSEVec[:,3]

bxGSEToCal    = bxGSEVec[timeFgmArrToCal]
byGSEToCal    = byGSEVec[timeFgmArrToCal]
bzGSEToCal    = bzGSEVec[timeFgmArrToCal]
btGSEToCal    = btGSEVec[timeFgmArrToCal]
epochFgmToCal = epochFgm[timeFgmArrToCal]
bGSEToCal     = np.column_stack((bxGSEToCal, byGSEToCal, bzGSEToCal))#, btGSEToCal))

[lambda1, lambda2, lambda3, N, M, L, eig, arrRot] = mva(bxGSEToCal, byGSEToCal, bzGSEToCal)
bGSEToRot = np.column_stack((bxGSEToCal, byGSEToCal, bzGSEToCal))
bLMN = np.dot(bGSEToRot, arrRot)

xlimEpochFgm = [min(epochFgmToCal), max(epochFgmToCal)]
colors    = ['blue', 'green', 'red']#, 'black']
labelsFgm = ['bx gse', 'by gse', 'bz gse']#, 'Bt GSE']
labelsLMN = ['bn', 'bm', 'bl']

kLoop     = np.arange(0,3)
for k in kLoop:
    ax0.plot(epochFgmToCal, bGSEToCal[:,k], linewidth = 0.8,
             color = colors[k], label=labelsFgm[k])
    ax0.legend(loc = 2, fontsize = 6, markerscale = 1.0,
               handlelength = 2, handleheight = 0.8)
'''
ax1.plot(epochFgmToCal, bLMN[:, 0], linewidth = 0.8,
         color = colors[0], label = labelsLMN[0])
'''
ax1.plot(epochFgmToCal, (-1)*bLMN[:, 0], linewidth = 0.8,
         color = colors[0], label = labelsLMN[0])
ax1.plot(epochFgmToCal, bLMN[:, 1], linewidth = 0.8,
         color = colors[1], label = labelsLMN[1])

ax1.plot(epochFgmToCal, (-1) * bLMN[:,2], linewidth = 0.8,
         color = colors[2], label = labelsLMN[2])
'''
ax1.plot(epochFgmToCal, bLMN[:,2], linewidth = 0.8,
         color = colors[2], label = labelsLMN[2])
'''
ax1.legend(loc = 2, fontsize = 6, markerscale = 0.5,
           handlelength = 2, handleheight = 0.8)
ax0.set_xlim(xlimEpochFgm)
ax1.set_xlim(xlimEpochFgm)
ax0.set_ylabel('nT')
ax1.set_ylabel('nT')

xtick0 = [2017, 12, 1, 14, 40, 50, 000]
xtick1 = [2017, 12, 1, 14, 40, 55, 000]
xtick2 = [2017, 12, 1, 14, 41,  0, 000]
xtick3 = [2017, 12, 1, 14, 41,  5, 000]

epochFgmTick = [xtick0, xtick1, xtick2, xtick3]
epochFgmTick = cdflib.cdfepoch.compute_tt2000(epochFgmTick)

###=================================================================================================
fileNameDisMom  = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-moms/2017/12/01/mms1_fpi_brst_l2_dis-moms_20171201143933_v3.3.0.cdf')
infoDisMom      = fileNameDisMom.cdf_info()
epochDisMom     = fileNameDisMom.varget('Epoch')
paraTVec        = fileNameDisMom.varget('mms1_dis_temppara_brst')
perpTVec        = fileNameDisMom.varget('mms1_dis_tempperp_brst')
enerSpecOmniVec = fileNameDisMom.varget('mms1_dis_energyspectr_omni_brst')
numDenVec       = fileNameDisMom.varget('mms1_dis_numberdensity_brst')

epochDisMomStartType2000 = cdflib.cdfepoch.compute_tt2000(timeInput0)
epochDisMomStopType2000  = cdflib.cdfepoch.compute_tt2000(timeInput1)
timeDisMomTuple0         = np.where(epochDisMom > epochDisMomStartType2000)
timeDisMomTuple1         = np.where(epochDisMom < epochDisMomStopType2000)
timeDisMomArr0           = timeDisMomTuple0[0]
timeDisMomArr1           = timeDisMomTuple1[0]
timeDisMomArrToCal       = np.intersect1d(timeDisMomArr0, timeDisMomArr1)

epochDisMomToCal  = epochDisMom[timeDisMomArrToCal]
paraTToCal        = paraTVec[timeDisMomArrToCal]
perpTToCal        = perpTVec[timeDisMomArrToCal]
tempToCal         = 1/3 * paraTToCal + 2/3 * perpTToCal###!!!!!!!!!!!!
enerSpecOmniToCal = enerSpecOmniVec[timeDisMomArrToCal]
numDenToCal       = numDenVec[timeDisMomArrToCal]

xlimEpochDisMom = [min(epochDisMomToCal), max(epochDisMomToCal)]
ax2.set_xlim(xlimEpochDisMom)
ax3.set_xlim(xlimEpochDisMom)
ax4.set_xlim(xlimEpochDisMom)
ax5.set_xlim(xlimEpochDisMom)
ax2.set_ylabel('$cm^{-1}$')

ax2.plot(epochDisMomToCal, numDenToCal,linewidth = 0.8,
         color = 'green', label = 'N ion')
ax2.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)
ax2.set_ylim([0,50])


ax3.plot(epochFgmToCal, btGSEToCal, linewidth = 0.8,
         color = 'black', label = 'Bt GSE')
ax3.legend(loc = 6, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)
ax3.set_ylim([0,40])
ax32 = ax3.twinx()
ax3.set_ylabel('nT')
ax32.plot(epochDisMomToCal, perpTToCal, linewidth = 0.8,
         color = 'steelblue', label = 'T perp')
ax32.legend(loc = 0, fontsize = 6, markerscale = 1.0,
           handlelength = 1, handleheight = 0.8)
ax32.set_ylabel('eV')
ax32.set_ylim([0,1600])

ax4.plot(epochDisMomToCal, paraTToCal, linewidth = 0.8,
         color = 'darkorange', label = 'T para')
ax4.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)
ax4.set_ylabel('eV')
ax4.set_ylim([0,1000])

ax5.plot(epochDisMomToCal, tempToCal,  linewidth = 0.8,
         color = 'black', label = 'T')
ax5.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)
ax5.set_ylabel('eV')
ax5.set_ylim([0,1500])

xtick0 = [2017, 12, 1, 14, 40, 50, 000]
xtick1 = [2017, 12, 1, 14, 40, 55, 000]
xtick2 = [2017, 12, 1, 14, 41,  0, 000]
xtick3 = [2017, 12, 1, 14, 41,  5, 000]

epochDisMomTick = [xtick0, xtick1, xtick2, xtick3]
epochDisMomTick = cdflib.cdfepoch.compute_tt2000(epochDisMomTick)
strDisMomTick = ['14:40:50', '14:40:55', '14:41:00', '14:41:05']
ax5.set_xticks(epochDisMomTick)
ax5.set_xticklabels(strDisMomTick, fontsize = 12)






pyEndTime = time.time()
print(
    '=================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))
###=====================================================================================================================
fileNameDesMom = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/des-moms/2017/12/01/mms1_fpi_brst_l2_des-moms_20171201143933_v3.3.0.cdf')
infoDesMom     = fileNameDesMom.cdf_info()
epochDesMom    = fileNameDesMom.varget('Epoch')
paraEleTVec    = fileNameDesMom.varget('mms1_des_temppara_brst')
perpEleTVec    = fileNameDesMom.varget('mms1_des_tempperp_brst')
numEleDenVec   = fileNameDesMom.varget('mms1_des_numberdensity_brst')

epochDesMomStartType2000 = cdflib.cdfepoch.compute_tt2000(timeInput0)
epochDesMomStopType2000  = cdflib.cdfepoch.compute_tt2000(timeInput1)
timeDesMomTuple0         = np.where(epochDesMom > epochDesMomStartType2000)
timeDesMomTuple1         = np.where(epochDesMom < epochDesMomStopType2000)
timeDesMomArr0           = timeDesMomTuple0[0]
timeDesMomArr1           = timeDesMomTuple1[0]
timeDesMomArrToCal       = np.intersect1d(timeDesMomArr0, timeDesMomArr1)

epochDesMomToCal = epochDesMom[timeDesMomArrToCal]
paraEleTToCal    = paraEleTVec[timeDesMomArrToCal]
perpEleTToCal    = perpEleTVec[timeDesMomArrToCal]
tempEleToCal     = 1/3*paraEleTToCal + 2/3*perpEleTToCal#####！！！！！！！！！！！！！！！！！！
numEleDenToCal   = numEleDenVec[timeDesMomArrToCal]
'''
ax2.plot(epochDesMomToCal, numEleDenToCal, linewidth = 0.8,
         color = 'lime', label = 'N ele')
ax2.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)

ax32.plot(epochDesMomToCal, perpEleTToCal, linewidth = 0.8,
          color = 'blue', label = 'T ele perp')
ax32.legend(loc = 0, fontsize = 6, markerscale = 1.0,
            handlelength = 1, handleheight = 0.8)

ax4.plot(epochDesMomToCal, paraEleTToCal, linewidth = 0.8,
         color = 'orange', label = 'T ele para')
ax4.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)

ax5.plot(epochDesMomToCal, tempEleToCal, linewidth = 0.8,
         color = 'grey', label = 'T ele')
ax5.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)
'''

plt.savefig(figurePath + '\\fluxRopeOrNot.png', dpi=199)
plt.show()
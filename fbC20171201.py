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

pyStartTime = time.time()

figurePath = 'D:\\wmmFigure\\FBs\\20171201\\'
if not os.path.exists(figurePath):
    os.makedirs(figurePath)



timeInput0 = [2017, 12, 1, 14, 40, 10, 000]
timeInput1 = [2017, 12, 1, 14, 41, 40, 000]

fig1 = plt.figure(figsize = (8, 10), dpi = 200)
fig1.subplots_adjust(top = 0.94)
fig1.subplots_adjust(bottom = 0.07)
fig1.subplots_adjust(left = 0.13)
fig1.subplots_adjust(right = 0.93)
gs = GridSpec(6,1, figure = fig1)
ax0 = fig1.add_subplot(gs[0,:])
ax1 = fig1.add_subplot(gs[1,:])
ax2 = fig1.add_subplot(gs[2,:])
ax3 = fig1.add_subplot(gs[3,:])
ax4 = fig1.add_subplot(gs[4,:])
ax5 = fig1.add_subplot(gs[5,:])

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
colors       = ['blue', 'green', 'red']#, 'black']
labelsFgm    = ['Bx GSE', 'By GSE', 'Bz GSE']#, 'Bt GSE']
kLoop        = np.arange(0,3)
for k in kLoop:
    ax0.plot(epochFgmToCal, bGSEToCal[:,k], linewidth = 0.8,
             color = colors[k], label=labelsFgm[k])
    ax0.legend(loc = 2, fontsize = 6, markerscale = 1.0,
               handlelength = 2, handleheight = 0.8)

xlimEpochFgm = [min(epochFgmToCal), max(epochFgmToCal)]

xtick0 = [2017, 12, 1, 14, 40, 15, 000]
xtick1 = [2017, 12, 1, 14, 40, 30, 000]
xtick2 = [2017, 12, 1, 14, 40, 45, 000]
xtick3 = [2017, 12, 1, 14, 41,  0, 000]
xtick4 = [2017, 12, 1, 14, 41, 15, 000]
xtick5 = [2017, 12, 1, 14, 41, 30, 000]

xtick0 = cdflib.cdfepoch.compute_tt2000(xtick0)
xtick1 = cdflib.cdfepoch.compute_tt2000(xtick1)
xtick2 = cdflib.cdfepoch.compute_tt2000(xtick2)
xtick3 = cdflib.cdfepoch.compute_tt2000(xtick3)
xtick4 = cdflib.cdfepoch.compute_tt2000(xtick4)
xtick5 = cdflib.cdfepoch.compute_tt2000(xtick5)

epochFgmTick = [xtick0, xtick1, xtick2, xtick3, xtick4, xtick5]
strFgmTick   = ['14:40:15', '14:40,30', '14:40:45', '14:41:00',
                '14:41:15', '14:41:30']

ax0.set_xlim(xlimEpochFgm)
ax0.set_ylim([-25, 35])
ax0.set_xticks(epochFgmTick)
ax0.set_xticklabels(strFgmTick, fontsize = 6)


ax1.plot(epochFgmToCal, btGSEToCal, linewidth = 0.8, color = 'black',
         label = 'Bt GSE')
ax1.legend(loc = 2, fontsize = 6, markerscale = 1.0, handlelength = 2, handleheight = 0.8)
ax1.set_xlim(xlimEpochFgm)
#ax1.set_ylim([])
ax1.set_xticks(epochFgmTick)
ax1.set_xticklabels(strFgmTick, fontsize = 6)
###=================================================================================================


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
tempToCal         = 2/3 * paraTToCal + 1/3 * perpTToCal
enerSpecOmniToCal = enerSpecOmniVec[timeDisMomArrToCal]
numDenToCal       = numDenVec[timeDisMomArrToCal]

xlimEpochDisMom = [min(epochDisMomToCal), max(epochDisMomToCal)]

ax2.plot(epochDisMomToCal, numDenToCal, linewidth = 0.8,
         color = 'green', label = 'N ion')
ax2.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)
ax2.set_xlim(xlimEpochDisMom)
ax2.set_xticks(epochFgmTick)
ax2.set_xticklabels(strFgmTick, fontsize = 6)


enerSpecOmniToCalT  = np.transpose(enerSpecOmniToCal)
#enerSpecOmniToCalTR = enerSpecOmniToCalT[::-1]
#im = ax5.imshow(enerSpecOmniToCalTR, aspect = 'auto', alpha=  0.8,
#                vmax = 0.7 * enerSpecOmniToCalTR.max(),
#                vmin = 0.3 * enerSpecOmniToCalTR.min())
imEnerSpec = ax3.imshow(enerSpecOmniToCalT,
#                       norm = colors.Normalize(vmin=-1.0, vmax=1.0),
                        vmin = 0.3 * enerSpecOmniToCalT.min(),
                        vmax = 0.7 * enerSpecOmniToCalT.max(),
                        alpha = 0.8,
                        aspect = 'auto',
                        origin = 'lower',
                        cmap = 'jet')
posCBar = [0.95,0.37, 0.01,0.12]
cbarEnerSpec = fig1.add_axes(posCBar)
plt.colorbar(imEnerSpec, cax=cbarEnerSpec)
cbarEnerSpec.set_label('energy flux')
ax3.set_yticks([])
ax3.set_xticks([])


ax4.plot(epochDisMomToCal, perpTToCal, linewidth = 0.8,
         color = 'steelblue', label = 'T perp')
ax4.plot(epochDisMomToCal, paraTToCal, linewidth = 0.8,
         color ='darkorange', label = 'T para')
ax4.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)

ax4.set_xlim(xlimEpochDisMom)
ax4.set_ylim([0, 2000])
ax4.set_xticks(epochFgmTick)
ax4.set_xticklabels(strFgmTick, fontsize = 6)

ax5.plot(epochDisMomToCal, tempToCal,  linewidth = 0.8,
         color = 'black',     label = 'T')
ax5.legend(loc = 2, fontsize = 6, markerscale = 1.0,
           handlelength = 2, handleheight = 0.8)
ax5.set_xlim(xlimEpochDisMom)
ax5.set_ylim([0, 2000])
ax5.set_xticks(epochFgmTick)
ax5.set_xticklabels(strFgmTick, fontsize = 6)








plt.savefig(figurePath + '\\fb.png', dpi=200)
plt.show()


###==================================================================
u,s,vh = np.linalg.svd(enerSpecOmniToCalT)

fig2 = plt.figure(figsize = (8, 10), dpi = 100)
fig2.subplots_adjust(top = 0.94)
fig2.subplots_adjust(bottom = 0.07)
fig2.subplots_adjust(left = 0.13)
fig2.subplots_adjust(right = 0.93)
gs = GridSpec(6,1, figure = fig2)
ax0 = fig2.add_subplot(gs[0,:])
ax1 = fig2.add_subplot(gs[1,:])
ax2 = fig2.add_subplot(gs[2,:])
ax3 = fig2.add_subplot(gs[3,:])
ax4 = fig2.add_subplot(gs[4,:])
ax5 = fig2.add_subplot(gs[5,:])

ax1.imshow(u)
ax2.imshow(vh)

plt.savefig(figurePath + '\\waveSVD.png', dpi=200)
plt.show()

pyEndTime = time.time()
print(
    '==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))

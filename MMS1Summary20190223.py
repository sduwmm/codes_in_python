#NormalOuter = [0.90525003, 0.076648012, 0.41790843]
#NormalInner=[0.91380150,-0.34255629,0.21822466]

# -*- coding: utf-8 -*-
'''
Created on Sat Apr 13 18:36:34 2019
Modified on Sun Apr 14
Changed on Mar 29 14:58:23 2021

@author: Mengmeng Wang
'''

import cdflib
import numpy as np
import math
import pyquaternion
import matplotlib.pyplot as plt
import os

time1=[2019,2,23,5,22,4,0]
time2=[2019,2,23,5,22,6,500]
time3=[2019,2,23,5,22,9,000]
time4=[2019,2,23,5,22,10,800]
time5=[2019,2,23,5,22,12,0]
time6=[2019,2,23,5,22,13,400]###Discontinuity
time7=[2019,2,23,5,22,14,500]###Discontinuity
time8=[2019,2,23,5,22,15,0]###Discontinuity
time9=[2019,2,23,5,22,18,0]
time10=[2019,2,23,5,22,20,0]

epoch1 = cdflib.cdfepoch.compute_tt2000(time1)
epoch2 = cdflib.cdfepoch.compute_tt2000(time2)
epoch3 = cdflib.cdfepoch.compute_tt2000(time3)
epoch4 = cdflib.cdfepoch.compute_tt2000(time4)
epoch5 = cdflib.cdfepoch.compute_tt2000(time5)
epoch6 = cdflib.cdfepoch.compute_tt2000(time6)
epoch7 = cdflib.cdfepoch.compute_tt2000(time7)
epoch8 = cdflib.cdfepoch.compute_tt2000(time8)
epoch9 = cdflib.cdfepoch.compute_tt2000(time9)
epoch10= cdflib.cdfepoch.compute_tt2000(time10)
epochList = [epoch1, epoch2, epoch3, epoch4, epoch5, epoch6, epoch7,epoch8, epoch9, epoch10]
epochArr  = np.array(epochList)


pthSW = [0.0334,0.0393,0.0358,0.0044,0.0140,0.0343,0,0.0676,0.0171,0.0172]
pdynSW=[0,0,0,0,0,0.4693,0,1.9687,0,0]
pthRE = [0.0907,0.0631,0.0846,0.2091,0.2504,0.1977,0,0.0814,0.0770,0.0843]
pdynRE=[0,0,0,0,0,0.0693,0,0.0052,0,0]
pb0=[0.0044,0.0056,0.0021,0.0015,0.0093,0.0421,0.0869,0.0425,0.0044,0.0048]
pthSWArr=np.array(pthSW)
pthREArr=np.array(pthRE)
pthAllArr=pthSWArr+pthREArr
ptotalSWArr=pthSWArr+np.array(pb0)+np.array(pdynSW)
ptotalREArr=pthREArr+np.array(pb0)+np.array(pdynRE)
ptotalAllArr=pthSWArr+np.array(pdynSW)+pthREArr+np.array(pdynRE)+np.array(pb0)
ptotalSWArr[6]=0
ptotalREArr[6]=0
pthAllArr[6]=0.2837
ptotalAllArr[6]=2.341940802132781




timeStartInput = [2019, 2, 23, 5, 22, 3, 000]
timeStopInput  = [2019, 2, 23, 5, 22, 22, 000]

epochStartType2000 = cdflib.cdfepoch.compute_tt2000(timeStartInput)
epochStopType2000  = cdflib.cdfepoch.compute_tt2000(timeStopInput)

fileNameFgm = cdflib.CDF('D:\\data\\mms\\mms1\\fgm\\brst\\l2\\2019\\02\\23\\mms1_fgm_brst_l2_20190223052123_v5.180.0.cdf')
infoFgm     = fileNameFgm.cdf_info()
epochFgm    = fileNameFgm.varget('Epoch')
bGSEVec     = fileNameFgm.varget('mms1_fgm_b_gse_brst_l2')

timeFgmTuple1    = np.where(epochFgm > epochStartType2000)
timeFgmTuple2    = np.where(epochFgm < epochStopType2000)
timeFgmArr1      = timeFgmTuple1[0]
timeFgmArr2      = timeFgmTuple2[0]
timeFgmArrToPlot = np.intersect1d(timeFgmArr1, timeFgmArr2)
epochFgmToPlot   = epochFgm[timeFgmArrToPlot]


fileNameMom = cdflib.CDF('D:\\data\\mms\\mms1\\fpi\\brst\\l2\\dis-moms\\2019\\02\\23\\mms1_fpi_brst_l2_dis-moms_20190223052123_v3.3.0.cdf')
infoMom     = fileNameMom.cdf_info()
epochMom    = fileNameMom.varget('Epoch')
densityMom  = fileNameMom.varget('mms1_dis_numberdensity_brst')


timeMomTuple1    = np.where(epochMom > epochStartType2000)  ###This is tuple data
timeMomTuple2    = np.where(epochMom < epochStopType2000)  ###This is tuple data
timeMomArr1      = timeMomTuple1[0]
timeMomArr2      = timeMomTuple2[0]
timeMomArrToPlot = np.intersect1d(timeMomArr1, timeMomArr2)
epochMomToPlot   = epochMom[timeMomArrToPlot]


bxArr  = bGSEVec[:,0]
byArr  = bGSEVec[:,1]
bzArr  = bGSEVec[:,2]
btArr  = bGSEVec[:,3]
denArr = densityMom


bxArrToPlot  = bxArr[timeFgmArrToPlot]
byArrToPlot  = byArr[timeFgmArrToPlot]
bzArrToPlot  = bzArr[timeFgmArrToPlot]
btArrToPlot  = btArr[timeFgmArrToPlot]
denArrToPlot = denArr[timeMomArrToPlot]

mu0  = 4*math.pi*10**(-7)
pb    = (btArrToPlot*10**(-9))**2/(2*mu0)
pbnPa = pb*10**9###magnetic pressure in nPa


fig = plt.figure(figsize=(16,20))#,dpi=500)
fig.subplots_adjust(top = 0.94)
fig.subplots_adjust(bottom = 0.07)
fig.subplots_adjust(left = 0.13)
fig.subplots_adjust(right = 0.93)
ax1 = plt.subplot(7,1,1)
ax2 = plt.subplot(7,1,2)
ax3 = plt.subplot(7,1,3)
ax4 = plt.subplot(7,1,4)
ax5 = plt.subplot(7,1,5)
ax6 = plt.subplot(7,1,6)
ax7 = plt.subplot(7,1,7)
plt.subplots_adjust(wspace =0, hspace =0.0)
#fig.tight_layout()


ax1.set_title('MMS 1 Observation: 2019-02-23', fontsize = 27)


xlimFgm = [min(epochFgmToPlot), max(epochFgmToPlot)]
xlimMom = [min(epochMomToPlot), max(epochMomToPlot)]
ax1.plot(epochFgmToPlot, bxArrToPlot, color = 'blue',  label = 'Bx GSE')
ax1.plot(epochFgmToPlot, byArrToPlot, color = 'green', label = 'By GSE')
ax1.plot(epochFgmToPlot, bzArrToPlot, color = 'red',   label = 'Bz GSE')
ax1.legend(loc = 2, fontsize = 15, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
#ax1.set_ylim([-15, 15])
ax1.set_xlim(xlimFgm)
ax1.set_xticks([])
ax1.tick_params(labelsize=15)
ax1.set_ylabel('B field [nT]', fontsize = 18)
ax1.text(epochFgmToPlot[2150], 2, '(a)', fontsize = 22)



ax2.plot(epochFgmToPlot, btArrToPlot,  color = 'black')
ax2.set_xlim(xlimFgm)
#ax2.set_ylim([0, 20])
ax2.set_xticks([])
ax2.set_ylabel('B total [nT]', color = 'black', fontsize = 18)
ax2.tick_params(axis = 'y', labelcolor = 'black', labelsize=15)
ax2.text(epochFgmToPlot[2150], 6, '(b)', fontsize = 22)

ax22 = ax2.twinx()
ax22.plot(epochMomToPlot, denArrToPlot, color = 'green')
#ax22.set_ylim([0, 20])
ax22.set_ylabel(r'$N [cm^{-3}]$', color = 'green', fontsize = 18)
ax22.tick_params(axis = 'y', labelcolor = 'green', labelsize = 15)


arrZero  = np.zeros((1, epochFgmToPlot.size))
arrZero  = arrZero.ravel()



#ax3.plot(epochMomToPlot, arrZero, 'k', linewidth = 0)
ax3.plot(epochFgmToPlot, pbnPa, color = 'black', label = 'pb', linewidth = 2)
ax3.scatter(epochArr, pthSWArr, color = 'orangered', label = 'pth sw', linewidth = 2)
ax3.legend(loc = 2, fontsize = 15, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax3.set_ylabel('P [nPa]', fontsize = 18)
ax3.set_ylim([0,0.4])
ax3.set_xlim(xlimFgm)
ax3.tick_params(labelsize=15)
ax3.set_xticks([])
ax3.text(epochFgmToPlot[2150], 0.12, '(c)', fontsize = 22)




ax4.plot(epochFgmToPlot, pbnPa,  color = 'black',   label = 'pb', linewidth = 2)
ax4.scatter(epochArr, pthREArr,  color = 'deepskyblue', label = 'pth re', linewidth = 2)
ax4.set_ylabel('P [nPa]', fontsize = 18)
ax4.legend(loc = 2, fontsize = 15, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax4.set_xlim(xlimMom)
ax4.set_ylim([0, 0.4])
ax4.tick_params(labelsize=15)
ax4.set_xticks([])
ax4.text(epochFgmToPlot[2150], 0.115, '(d)', fontsize = 22)



ax5.plot(epochFgmToPlot, arrZero, 'k', linewidth = 0)
ax5.scatter(epochArr, pthAllArr,  color = 'crimson', marker = 'o', linewidth = 2, label = 'pth all')
ax5.plot(epochFgmToPlot, pbnPa, color = 'black', label = 'pb', linewidth = 2)
ax5.legend(loc = 2, fontsize = 15, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax5.set_ylabel('p [nPa]', fontsize = 18)
ax5.set_xlim(xlimFgm)
ax5.set_ylim([0, 0.4])
ax5.tick_params(labelsize=15)
ax5.set_xticks([])
ax5.text(epochFgmToPlot[2150], 0.155, '(e)', fontsize = 22)



ax6.plot(epochFgmToPlot, arrZero, 'k', linewidth = 0)
ax6.scatter(epochArr, ptotalSWArr, color = 'orangered',  marker = '*', label = 'Solar wind beam', linewidth = 3)
ax6.scatter(epochArr, ptotalREArr, color = 'deepskyblue', marker = '*', label = 'Reflected ion beam', linewidth = 3)
ax6.legend(loc = 2, fontsize = 15, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax6.set_ylabel('pt [nPa]', fontsize = 18)
ax6.set_xlim(xlimFgm)
ax6.set_ylim([0, 3])
ax6.tick_params(labelsize=15)
ax6.set_xticks([])
ax6.text(epochFgmToPlot[2150], 0.155, '(f)', fontsize = 22)



xtick1 = [2019, 2, 23, 5, 22, 5,0]
xtick2 = [2019, 2, 23, 5, 22, 10,0]
xtick3 = [2019, 2, 23, 5, 22, 15,0]
xtick4 = [2019, 2, 23, 5, 22, 20,0]

xtick1 = cdflib.cdfepoch.compute_tt2000(xtick1)
xtick2 = cdflib.cdfepoch.compute_tt2000(xtick2)
xtick3 = cdflib.cdfepoch.compute_tt2000(xtick3)
xtick4 = cdflib.cdfepoch.compute_tt2000(xtick4)

epochTick = [xtick1, xtick2, xtick3, xtick4 ]
strTick = ['5:22:5', '5:22:10', '5:22:15', '5:22:20']




ax7.plot(epochFgmToPlot, arrZero, 'k', linewidth = 0)
ax7.scatter(epochArr, ptotalAllArr, color = 'black', marker = 'D', linewidth = 2, label = 'all')
ax7.legend(loc = 2, fontsize = 15, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax7.set_ylabel('pt [nPa]', fontsize = 18)
ax7.set_xlim(xlimFgm)
ax7.set_ylim([0, 3])
ax7.tick_params(labelsize=15)
ax7.set_xticks([])
ax7.text(epochFgmToPlot[2150], 0.24, '(g)', fontsize = 22)
ax7.set_xticks(epochTick)
ax7.set_xticklabels(strTick, fontsize = 25)
ax7.set_xticks(epochTick, strTick)


figurePath = 'D:/OneDrive - mail.sdu.edu.cn/work/code&figure/figure2021/'
if not os.path.exists(figurePath):
    os.makedirs(figurePath)

plt.savefig(figurePath+'/MMS1Pressure.png',dpi = 200)
plt.show()
#plt.close()
#pdynSWList = [2.0945, 0.7037, 0.2977, 1.3954, 1.6708, 0.8505, 2.5500]

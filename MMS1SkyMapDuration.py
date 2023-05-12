
# -*- coding: utf-8 -*-
"""
Created on Thur Mar 18 15:21:30 2021
@author: Mengmeng Wang
"""
"""
plot the sky map
of MMS FPI ion distribution data
from the phase space density

Version 0
"""

"""
!!!Note
dist: time*energy*theta*phi=epochNumber*32*16*32
"""

import os
import cdflib
import matplotlib.pyplot as plt
import numpy as np
import math
import time

pyStartTime = time.time()

TimeStart = [2018, 3, 4, 16, 30, 0, 0]
TimeEnd = [2018, 3, 4, 16, 30, 5, 0]
FileDist = \
    cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2018/03/04/mms1_fpi_brst_l2_dis-dist_20180304162923_v3.3.0.cdf')

FigurePath = 'D:/OneDrive - mail.sdu.edu.cn/work/code&figure/figure2021/skyMap/'
#
InfoDist = FileDist.cdf_info()
Epoch = FileDist.varget('Epoch')
Phi = FileDist.varget('mms1_dis_phi_brst')
PhiDelta = FileDist.varget('mms1_dis_theta_brst')
Dist = FileDist.varget('mms1_dis_dist_brst')
DistErr = FileDist.varget('mms1_dis_disterr_brst')
Theta = FileDist.varget('mms1_dis_theta_brst')
ThetaDelta = FileDist.varget('mms1_dis_theta_delta_brst')
Energy = FileDist.varget('mms1_dis_energy_brst')
EnergyDelta = FileDist.varget('mms1_dis_energy_delta_brst')
#
EpochStartType2000 = cdflib.cdfepoch.compute_tt2000(TimeStart)
EpochEndType2000 = cdflib.cdfepoch.compute_tt2000(TimeEnd)
TimeTplStart = np.where(Epoch >= EpochStartType2000)  ###This is tuple data
TimeTplEnd = np.where(Epoch <= EpochEndType2000)  ###This is tuple data
TimeArr1 = TimeTplStart[0]
TimeArr2 = TimeTplEnd[0]
TimeArrToCal = np.intersect1d(TimeArr1, TimeArr2)
#
for tt in TimeArrToCal:
    EnergyTT = Energy[tt]
    DistTT = Dist[tt]
    DistErrTT = Dist[tt]
    DistAlongEnergySum = np.zeros([16, 32])
    for ee in np.arange(32):
        EnergyInd = ee
        DistTTEE = DistTT[ee]#take out the dist of certain energy channel
        DistErrTTEE = DistErrTT[ee]#take out the dist error of certain energy channel
        #
        #DistTTEE[DistTTEE <= DistErrTTEE] = 0
        #Ind = np.where(DistTTEE <= DistErrTTEE)
        #DistTTEE[Ind] = 0
        #
        DistAlongEnergySum = DistAlongEnergySum + DistTTEE
        #
        DistSumLg = np.log10(DistAlongEnergySum)
    #
    Figure = plt.figure(figsize = (5,4), dpi = 100)
    ax = plt.subplot(1,1,1)
    cmap = ax.pcolor(DistSumLg, cmap = 'jet', vmin= -26, vmax = -20)
    cbar = Figure.colorbar(cmap, ax = ax)
    cbar.set_label(r'$Log_{10} PSD (s^3/cm^6)$', x=0, y=0.55, rotation=90, fontsize=10)
    ax.text(4, 16.5, '2018-02-23/16:30:00+' + str(tt-TimeArrToCal[0]), fontsize=10)
    #
    if not os.path.exists(FigurePath):
        os.makedirs(FigurePath)
    FigureName = FigurePath + 'MMS1_' + str(tt-TimeArrToCal[0])

    plt.savefig(FigureName)
    plt.show()
    plt.close()
pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))


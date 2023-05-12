# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:13:30 2021
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


#TimeInput = [2018, 3, 4, 16, 30, 2, 312]
#FileName = \
#    cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2018/03/04/mms1_fpi_brst_l2_dis-dist_20180304162923_v3.3.0.cdf')

TimeInput = [2019, 2, 23, 5, 22, 5, 921]
FileName = \
    cdflib.CDF('D:/data/mms/mms2/fpi/brst/l2/dis-dist/2019/02/23/mms2_fpi_brst_l2_dis-dist_20190223052123_v3.3.0.cdf')
###This can become more automatable.
Info = FileName.cdf_info()
Epoch = FileName.varget('Epoch')
EpochPlusVar = FileName.varget('Epoch_plus_var')
Phi = FileName.varget('mms2_dis_phi_brst')
PhiDelta = FileName.varget('mms2_dis_phi_delta_brst')
Theta = FileName.varget('mms2_dis_theta_brst')
ThetaDelta = FileName.varget('mms2_dis_theta_delta_brst')
Dist = FileName.varget('mms2_dis_dist_brst')
DistErr = FileName.varget('mms2_dis_disterr_brst')
EnergyLabel = FileName.varget('mms2_dis_energy_label_brst')
Energy = FileName.varget('mms2_dis_energy_brst')
EnergyDelta = FileName.varget('mms2_dis_energy_delta_brst')

EpochToCalType2000 = cdflib.cdfepoch.compute_tt2000(TimeInput)
TimeInTpl = np.where(Epoch >= EpochToCalType2000)
TimeInArr = TimeInTpl[0]
TimeIndex = TimeInArr[0]

DistTimeEnerSum = np.zeros([16,32])
#for ee in np.arange(32):
for ee in np.arange(0,32,1):
    #
    EnergyIndex = ee
    #
    PhiTime = Phi[TimeIndex]
    ThetaTime = Theta
    DistTimeEner = Dist[TimeIndex, EnergyIndex]
    #
    DistTimeEnerSum = DistTimeEnerSum + DistTimeEner
    #
    DistTimeEnerSumLg = np.log10(DistTimeEnerSum)
EnergyPlotStart = Energy[TimeIndex, 0]
EnergyDeltaPlotStart = EnergyDelta[TimeIndex, 0]
EnergyPlotEnd   = Energy[TimeIndex,32-1]
EnergyDeltaPlotEnd = EnergyDelta[TimeIndex, 32-1]

Fig = plt.figure(figsize=(5, 4), dpi=200)
ax = plt.subplot(1, 1, 1)
#c = ax.pcolor(DistTimeEnerSumLg, cmap='tab20b', vmin=-26, vmax=-20)
c = ax.pcolor(DistTimeEnerSumLg, cmap='jet', vmin=-26, vmax=-20)
cbar = Fig.colorbar(c, ax=ax)
cbar.set_label(r'$Log_{10} PSD (s^3/cm^6)$', x=0, y=0.55, rotation=90, fontsize=10)

TimeStr = str(TimeInput[0]) + '-' + str(TimeInput[1]) + '-' + str(TimeInput[2]) + '/' + \
          str(TimeInput[3]) + ':' + str(TimeInput[4]) + ':' + str(TimeInput[5]) + '.' + \
          str(TimeInput[6])
# print(TimeStr)
print(EnergyPlotStart)
print(EnergyPlotEnd)
ax.text(4, 16.5, TimeStr, fontsize=10)
EnergyRangeStr = '{0:.0f}'.format(EnergyPlotStart - EnergyDeltaPlotStart) + '-' + \
                 '{0:.0f}'.format(EnergyPlotEnd + EnergyDeltaPlotEnd) + 'eV'
ax.text(21, 16.5, EnergyRangeStr, fontsize=10)

FigurePath = 'D:/OneDrive - mail.sdu.edu.cn/work/code&figure/figure2021/skyMap/'
if not os.path.exists(FigurePath):
    os.makedirs(FigurePath)

FigureName = FigurePath + str(TimeInput[0]) + str(TimeInput[1]) + str(TimeInput[2]) + \
             '-' + str(TimeInput[3]) + str(TimeInput[4]) + str(TimeInput[5]) + '.' + str(TimeInput[6]) + \
             '.png'

plt.savefig(FigureName)
plt.show()
plt.close()
# 'Energy' + str(EnergyIndex) + '.png'  # {:-5}.png', format(EnergyIndex)



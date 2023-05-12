# -*- coding: utf-8 -*-
"""
Created on Thur Apr 4 17:04:37 2019
Modified on Tue Apr 9 11:03:04 2019
@author: Mengmeng Wang
"""

"""
plot the sky map
of MMS FPI ion distribution data
from the phase space density

Verion 0
"""

"""
!!!Note
dist: time*energy*theta*phi=epoch_number*32*16*32
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import cdflib
import numpy as np
import math
import os

timeInput  = [2018, 3, 21, 21, 6, 14, 000]


fileName     = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2017/11/16/mms1_fpi_brst_l2_dis-dist_20171116121223_v3.3.0.cdf')
info         = fileName.cdf_info()
epoch        = fileName.varget('epoch')
epochPlusVar = fileName.varget('Epoch_plus_var')
phi          = fileName.varget('mms1_dis_phi_brst')
phiDelta     = fileName.varget('mms1_dis_phi_delta_brst')
theta        = fileName.varget('mms1_dis_theta_brst')
thetaDelta   = fileName.varget('mms1_dis_theta_delta_brst')
dist         = fileName.varget('mms1_dis_dist_brst')
distErr      = fileName.varget('mms1_dis_disterr_brst')
energyLabel  = fileName.varget('mms1_dis_energy_label_brst')
energy       = fileName.varget('mms1_dis_energy_brst')
energyDelta  = fileName.varget('mms1_dis_energy_delta_brst')

epochToCalType2000 = cdflib.cdfepoch.compute_tt2000(timeInput)
timeInTpl          = np.where(epoch > epochToCalType2000)  ###This is tuple data
timeInArr          = timeInTpl[0]
timeIndex          = timeInArr[0]
#timeToPlot = epoch[timeIndex]
for ee in np.arange(32):

    energyIndex     = ee
    energyPlot      = energy[timeIndex, ee]
    energyDeltaPlot = energyDelta[timeIndex, ee]

    phiTime      = phi[timeIndex]
    thetaTime    = theta
    distTimeEner = dist[timeIndex, energyIndex]
    #
    nanIndex               = np.where(distTimeEner == 0)
    distTimeEner[nanIndex] = np.NaN
    #
    distTimeEnerLg = np.log10(distTimeEner)

    #distTimeEner23 = distTimeEner*10**19

    #fig = plt.figure(0)
    fig  = plt.figure(figsize=(10, 8), dpi=200)
    ax   = plt.subplot(1, 1, 1)
    #fig, ax = plt.subplots(1,1)
    c    = ax.pcolor(distTimeEnerLg,  cmap='rainbow')
    cbar = fig.colorbar(c, ax = ax)
    cbar.set_label(r'$Log_{10} PSD(s^3/cm^6)$', x =0, y = 0.55, rotation=90, fontsize = 15)

    timeStr = str(timeInput[0]) + '-' + str(timeInput[1]) + '-' + str(timeInput[2]) + '/' + \
        str(timeInput[3]) + ':' + str(timeInput[4]) + ':' + str(timeInput[5]) + '.' + str(timeInput[6])
    ax.text(4, 16.5, timeStr, fontsize =15)
    energyRangeStr = '{0:.0f}'.format(energyPlot-energyDeltaPlot) + '-'+'{0:.0f}'.format(energyPlot+energyDeltaPlot) + 'eV'
    ax.text(21, 16.5, energyRangeStr, fontsize = 15)

    #ax.set_xticks([])
    #phiInt    = phiTime.astype(np.int32)
    #phiLabels = phiInt.tolist()
    #ax.set_xticks(phiTime)
    #ax.text(0, -1, phiLabels)

    #RdBu_r
    #RdYlBu_r
    #Spectral_r
    figurePath = 'D:/shfaFigure/20171116121300/skyMapAna/'
    if not os.path.exists(figurePath):
        os.makedirs(figurePath)

    figureName = figurePath + str(timeInput[0]) + str(timeInput[1]) + str(timeInput[2]) + '_' + str(timeInput[3]) +\
        str(timeInput[4]) +str(timeInput[5]) + '.'+str(timeInput[6])+\
        'energy' + str(energyIndex) + '.png'#{:-5}.png',format(energyIndex)
    plt.savefig(figureName)
    plt.show()
    plt.close()


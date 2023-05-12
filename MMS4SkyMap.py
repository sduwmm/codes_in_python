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

def PlotSkyMap(Time, FileDist, EnergyChannel, FigurePath):

    InfoDist = FileDist.cdf_info()
    Epoch = FileDist.varget('Epoch')
    Phi = FileDist.varget('mms4_dis_phi_brst')
    PhiDelta = FileDist.varget('mms4_dis_theta_brst')
    Dist = FileDist.varget('mms4_dis_dist_brst')
    DistErr = FileDist.varget('mms4_dis_disterr_brst')
    Theta = FileDist.varget('mms4_dis_theta_brst')
    ThetaDelta  = FileDist.varget('mms4_dis_theta_delta_brst')
    Energy = FileDist.varget('mms4_dis_energy_brst')
    EnergyDelta = FileDist.varget('mms4_dis_energy_delta_brst')

    EpochType2000 = cdflib.cdfepoch.compute_tt2000(Time)
    TimeIndTpl = np.where(Epoch >= EpochType2000)
    TimeIndArr = TimeIndTpl[0]
    TimeInd = TimeIndArr[0]
    ####initialize==========================
    DistAlongEnergySum = np.zeros([16,32])
    EnergyPlot = np.array(EnergyChannel)
    Theta = Theta
    PhiPlot = Phi[TimeInd]
    DistPlot = Dist[TimeInd]
    DistErrPlot = DistErr[TimeInd]
    ###==========================initialize
    for ee in EnergyPlot:
        EnergyInd = ee
        DistPlotEE = DistPlot[ee]#take out the dist of certain energy channel
        DistErrPlotEE = DistErrPlot[ee]#take out the dist error of certain energy channel
        #
        DistPlotEE[DistPlotEE <= DistErrPlotEE] = 0
        #
        DistAlongEnergySum =  DistAlongEnergySum + DistPlotEE
        #
        DistSumLg = np.log10(DistAlongEnergySum)
    EnergyStart = Energy[TimeInd, EnergyPlot[0]]
    EnergyDeltaStart = EnergyDelta[TimeInd, EnergyPlot[0]]
    EnergyEnd = Energy[TimeInd, EnergyPlot[-1]]
    EnergyDeltaEnd = Energy[TimeInd, EnergyPlot[-1]]
    #
    Figure = plt.figure(figsize = (5,4), dpi =100)
    ax = plt.subplot(1,1,1)
    cmap = ax.pcolor(DistSumLg, cmap = 'jet', vmin = -26, vmax = -20)
    cbar = Figure.colorbar(cmap, ax = ax)
    cbar.set_label(r'$Log_{10} PSD (s^3/cm^6)$', x=0, y=0.55, rotation=90, fontsize=10)
    #
    TimeStr = str(Time[0]) + '-' + str(Time[1]) + '-' + str(Time[2]) + '/' + \
              str(Time[3]) + ':' + str(Time[4]) + ':' + str(Time[5]) + '.' + \
              str(Time[6])
    # print(TimeStr)
    print(EnergyStart-EnergyDeltaStart)
    print(EnergyEnd+EnergyDeltaEnd)
    ax.text(4, 16.5, TimeStr, fontsize=10)
    EnergyRangeStr = '{0:.0f}'.format(EnergyStart - EnergyDeltaStart) + '-' + \
                     '{0:.0f}'.format(EnergyEnd + EnergyDeltaEnd) + 'eV'
    ax.text(21, 16.5, EnergyRangeStr, fontsize=10)

    if not os.path.exists(FigurePath):
        os.makedirs(FigurePath)
    FigureName = FigurePath + 'MMS4_' + str(Time[0]) + str(Time[1]) + str(Time[2]) + \
                 '-' + str(Time[3]) + str(Time[4]) + str(Time[5]) + '.' + str(Time[6]) + \
                 '.png'

    plt.savefig(FigureName)
    plt.show()
    #plt.close()
    return()
###========================================
TimeInput = [2019, 2, 23, 5, 22, 9, 921]
FileInput = \
    cdflib.CDF('D:/data/mms/mms4/fpi/brst/l2/dis-dist/2019/02/23/mms4_fpi_brst_l2_dis-dist_20190223052123_v3.3.0.cdf')

FigurePathInput = 'D:/OneDrive - mail.sdu.edu.cn/work/code&figure/figure2021/skyMap/'
EnergyInput = np.arange(0,32,1)
SkyMap = PlotSkyMap(TimeInput, FileInput, EnergyInput, FigurePathInput)
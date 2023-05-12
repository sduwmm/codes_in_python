# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 20:19:04 2019

@author: Mengmeng Wang
"""

"""
calculate the zero-, first-, second- moments
of MMS FPI ion distribution data
from the phase space density

Version test 
For reading cdf file
"""
"""
!!Note
dist: time*energy*theta*phi
"""

import cdflib
import numpy as np
import matplotlib.pyplot as plt
import datetime

fileName = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2017/11/16/mms1_fpi_brst_l2_dis-dist_20171116121223_v3.3.0.cdf')
info = fileName.cdf_info()
epoch = fileName.varget('epoch')
epochPlusVar = fileName.varget('Epoch_plus_var')
phi = fileName.varget('mms1_dis_phi_brst')
phiDelta = fileName.varget('mms1_dis_phi_delta_brst')
theta = fileName.varget('mms1_dis_theta_brst')
thetaDelta = fileName.varget('mms1_dis_theta_delta_brst')
dist = fileName.varget('mms1_dis_dist_brst')
distErr = fileName.varget('mms1_dis_disterr_brst')
energyLabel = fileName.varget('mms1_dis_energy_label_brst')
energy = fileName.varget('mms1_dis_energy_brst')
energyDelta = fileName.varget('mms1_dis_energy_delta_brst')
print(epoch)
print(type(epoch))
print(epoch.dtype)
print(epoch.size)
print(type(phi))
print(phi.dtype)
print(type(dist))
print(dist.dtype)
print(dist.ndim)
print(dist.size)
print(dist.shape)

print('Second part====================================')
timeToPlotInput = [2017,11,16,12,13,00,000]
epochToPlotType2000 = cdflib.cdfepoch.compute_tt2000(timeToPlotInput)
timeInTuple = np.where(epoch > epochToPlotType2000)##tuple data
timeInArr = timeInTuple[0]
timeIndex = timeInArr[0]
timeToPlot = epoch[timeIndex]

distFunTime = dist[timeIndex]
enerTime = energy[timeIndex]
print(distFunTime.shape)
print(distFunTime.data)
enerChannel = 31
distFunTimeEner = distFunTime[enerChannel,:,:]###time*energy*theta*phi
print(enerTime[enerChannel])
print(distFunTimeEner)
print(distFunTimeEner.size)

oneD = distFunTimeEner.ravel()
print(oneD)
print("========================================================")
oneDsort = np.sort(oneD, axis = 0)
print(oneDsort)
print(enerTime[enerChannel])
"""
print(epochToPlotType2000)
print(type(timeInTuple))
print(type(timeInArr))
print(timeToPlot)
"""
"""
print(distFunTime)
"""
# =============================================================================
# ctrl 4
#a = cdflib.cdfepoch.compute_tt2000([2017,11,16,12,13,30,150])
#b = cdflib.cdfepoch.compute_tt2000([2017,11,16,12,13,30,300])
#tt = b - a
#print(tt)
#A = [a,b,timeToPlotNew]
#B = np.sort(A, axis = 0)
#print(B)
# =============================================================================

"""
epoch[-1] - epoch[-2] = 150000000
"""
"""
print(info)
print(epoch)
print(epochPlusVar)
print(phi)
print(phiDelta)
print(theta)
print(thetaDelta)
print(dist)
print(distErr)
print(energyLabel)
print(energy)
print(energyDelta)
"""

"""
'Epoch'
'Epoch_plus_var'

'mms1_dis_phi_brst'
'mms1_dis_phi_delta_brst'

'mms1_dis_dist_brst'
'mms1_dis_disterr_brst'

'mms1_dis_theta_brst'
'mms1_dis_theta_delta_brst'

'mms1_dis_energy_label_brst'
'mms1_dis_energy_brst'
'mms1_dis_energy_delta_brst'
"""


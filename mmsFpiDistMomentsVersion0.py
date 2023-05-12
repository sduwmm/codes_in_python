# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 10:57:37 2019

@author: Mengmeng Wang
"""

"""
calculate the zero-, first-, second- moments
of MMS FPI ion distribution data
from the phase space density

Verion 0
"""

"""
!!!Note
dist: time*energy*theta*phi=epoch_number*32*16*32
"""

import cdflib
import numpy as np
import math
import sympy
# import scipy.integrate
import matplotlib.pyplot as plt

###constants==================
mp = 1.67262158 * 10 ** (-27)
eV2J = 1.602176462 * 10 ** (-19)
kb = 1.3806503 * 10 ** (-23)

timeInput = [2017, 11, 16, 12, 13, 00, 000]

fileName = cdflib.CDF(
    'D:/data/mms/mms1/fpi/brst/l2/dis-dist/2017/11/16/mms1_fpi_brst_l2_dis-dist_20171116121223_v3.3.0.cdf')
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

epochToCalType2000 = cdflib.cdfepoch.compute_tt2000(timeInput)
timeInTpl = np.where(epoch > epochToCalType2000)  ###This is tuple data
timeInArr = timeInTpl[0]
timeIndex = timeInArr[0]
timeToCal = epoch[timeIndex]

phiTime = phi[timeIndex] * math.pi / 180
thetaTime = theta * math.pi / 180
deltaPhi = phiDelta[0] * 2 * math.pi / 180
deltaTheta = thetaDelta[0] * 2 * math.pi / 180
distTime0 = dist[timeIndex]
distTime = distTime0 * 10 ** 12
enerTime = energy[timeIndex]
enerDelta = energyDelta[timeIndex]
###
enerTimeJ = enerTime * eV2J
enerDeltaJ = enerDelta * eV2J
###

loopEner = np.arange(0, enerTime.size)
loopPhi = np.arange(0, phiTime.size)
loopTheta = np.arange(0, thetaTime.size)

den = 0
for ii in loopEner:
    enerII = enerTimeJ[ii]
    enerDeltaII = enerDeltaJ[ii]
    velIISqr = 2 * enerII / mp
    enerIILeft = enerII - enerDeltaII
    enerIIRight = enerII + enerDeltaII
    velIILeft = math.sqrt(2 * enerIILeft / mp)
    velIIRight = math.sqrt(2 * enerIIRight / mp)
    velDelII = velIIRight - velIILeft
    for jj in loopTheta:
        thetaJJ = thetaTime[jj]
        for kk in loopPhi:
            denPlus = deltaTheta * deltaPhi * velIISqr * velDelII * math.sin(thetaJJ) * distTime[ii, jj, kk]
            den = den + denPlus

densityM = den  ###number density in m**(-3)
densityCM = den * 10 ** (-6)  ###number density in cm**(-3)

nByVx = 0
nByVy = 0
nByVz = 0
for ii in loopEner:
    enerII = enerTimeJ[ii]
    enerDeltaII = enerDeltaJ[ii]
    velII = math.sqrt(2 * enerII / mp)
    velIICube = velII ** 3
    enerIILeft = enerII - enerDeltaII
    enerIIRight = enerII + enerDeltaII
    velIILeft = math.sqrt(2 * enerIILeft / mp)
    velIIRight = math.sqrt(2 * enerIIRight / mp)
    velDelII = velIIRight - velIILeft
    for jj in loopTheta:
        thetaJJ = thetaTime[jj]
        for kk in loopPhi:
            phiKK = phiTime[kk]
            nByVxPlus = deltaTheta * deltaPhi * velIICube * velDelII * (math.sin(thetaJJ)) ** 2 * math.cos(phiKK) * \
                        distTime[ii, jj, kk]
            nByVx = nByVx + nByVxPlus
            nByVyPlus = deltaTheta * deltaPhi * velIICube * velDelII * (math.sin(thetaJJ)) * (math.sin(phiKK)) ** 2 * \
                        distTime[ii, jj, kk]
            nByVy = nByVy + nByVyPlus
            nByVzPlus = deltaTheta * deltaPhi * velIICube * velDelII * math.sin(thetaJJ) * math.cos(thetaJJ) * distTime[
                ii, jj, kk]
            nByVz = nByVz + nByVzPlus
vxMS = nByVx / densityM
vxKmS = vxMS * 10 ** (-3)
vyMS = nByVy / densityM
vyKmS = vyMS * 10 ** (-3)
vzMS = nByVz / densityM
vzKmS = vzMS * 10 ** (-3)

pxxItem1 = 0
pyyItem1 = 0
pzzItem1 = 0
for ii in loopEner:
    enerII = enerTimeJ[ii]
    enerDeltaII = enerDeltaJ[ii]
    velIISqr = 2 * enerII / mp
    velIIQuar = velIISqr ** 2
    enerIILeft = enerII - enerDeltaII
    enerIIRight = enerII + enerDeltaII
    velIILeft = math.sqrt(2 * enerIILeft / mp)
    velIIRight = math.sqrt(2 * enerIIRight / mp)
    velDelII = velIIRight - velIILeft
    for jj in loopTheta:
        thetaJJ = thetaTime[jj]
        for kk in loopPhi:
            phiKK = phiTime[kk]
            pxxPlus = mp * deltaTheta * deltaPhi * velIIQuar * velDelII * (math.sin(thetaJJ)) ** 3 * (
                math.cos(phiKK)) ** 2 * distTime[ii, jj, kk]
            pxxItem1 = pxxItem1 + pxxPlus
            pyyPlus = mp * deltaTheta * deltaPhi * velIIQuar * velDelII * (math.sin(thetaJJ)) ** 3 * (
                math.sin(phiKK)) ** 2 * distTime[ii, jj, kk]
            pyyItem1 = pyyItem1 + pyyPlus
            pzzPlus = mp * deltaTheta * deltaPhi * velIIQuar * velDelII * math.sin(thetaJJ) * (math.cos(thetaJJ)) ** 2 * \
                      distTime[ii, jj, kk]
            pzzItem1 = pzzItem1 + pzzPlus
pxx = pxxItem1 - mp * densityM * vxMS ** 2
pyy = pyyItem1 - mp * densityM * vxMS ** 2
pzz = pzzItem1 - mp * densityM * vzMS ** 2
p = (pxx + pyy + pzz) / 3
pInnPa = p * 10 ** 9
tempK = p / (densityM * kb)
tempEV = tempK / eV2J
print('density:', densityCM)
print('vx:', vxKmS)
print('vy:', vyKmS)
print('vz:', vzKmS)
print('pth nPa:', pInnPa)
print('temp eV:', tempEV)

'''
enerTimeNew = np.zeros(enerTime.size)
loopArr = np.arange(0,enerTime.size-1)
#for i in range(enerTime.size):
enerTimeNew[0] = enerTime[0]
for i in loopArr:
    enerTimeNew[i+1] = enerTime[i] + enerDelta[i] + enerDelta[i+1]

A = enerTimeNew - enerTime
plotArr = np.arange(0, enerTime.size)
fig, axs = plt.subplots(1,4,figsize = (12,3))
axs[0].bar(plotArr, enerTime)
axs[1].bar(plotArr, enerTimeNew)
axs[2].bar(plotArr, A)
'''

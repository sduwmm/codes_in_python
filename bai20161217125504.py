# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 09:48:28 2019
Modified on Thur Apr 9 16:31:40 2019

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
import matplotlib.pyplot as plt
###
import time


pyStartTime = time.time()
###
###constants==================
mp   = 1.67262158*10**(-27)
eV2J = 1.602176462*10**(-19)
kb   = 1.3806503*10**(-23)


timeInput1 = [2016,12,17,12,55, 4,000]
timeInput2 = [2016,12,17,12,55,43,000]

loopEnerO = np.arange(0,19)

fileName     = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2016/12/17/mms1_fpi_brst_l2_dis-dist_20161217125504_v3.3.0.cdf')
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

epochStartType2000 = cdflib.cdfepoch.compute_tt2000(timeInput1)
epochStopType2000  = cdflib.cdfepoch.compute_tt2000(timeInput2)
timeTuple1         = np.where(epoch > epochStartType2000)###This is tuple data
timeTuple2         = np.where(epoch < epochStopType2000)###This is tuple data
timeArr1           = timeTuple1[0]
timeArr2           = timeTuple2[0]
timeArrToCal       = np.intersect1d(timeArr1,timeArr2)

###convert data to new data in SI Units
phi        = phi*math.pi/180###angle in radian
theta      = theta*math.pi/180###angle in radian
deltaPhi   = phiDelta[0]*2*math.pi/180###angle in radian
deltaTheta = thetaDelta[0]*2*math.pi/180###angle in radian
dist0      = dist*10**12###dist in s3/m6
energyJ    = energy*eV2J###energy in J
enerDeltaJ = energyDelta*eV2J###energy in J
###convert data to new data in SI Units


den    = []
vx     = []
vy     = []
vz     = []
vtotal = []
pth    = []
for tt in timeArrToCal:
    phiTT       = phi[tt]
    thetaTT     = theta
    distTT      = dist0[tt]
    enerTT      = energyJ[tt]
    enerDeltaTT = enerDeltaJ[tt]

    #loopEner  = np.arange(0, enerTT.size)
    loopEner = loopEnerO
    loopPhi   = np.arange(0, phiTT.size)
    loopTheta = np.arange(0, thetaTT.size)

    denTT    = 0
    nByVxTT  = 0
    nByVyTT  = 0
    nByVzTT  = 0
    pxxItem1 = 0
    pyyItem1 = 0
    pzzItem1 = 0
    for ii in loopEner:
        enerII      = enerTT[ii]
        enerDeltaII = enerDeltaTT[ii]
        velII       = math.sqrt(2*enerII/mp)
        velIISqr    = velII**2
        velIICube   = velII**3
        velIIQuar   = velII**4
        enerIILeft  = enerII - enerDeltaII
        enerIIRight = enerII + enerDeltaII
        velIILeft   = math.sqrt(2*enerIILeft/mp)
        velIIRight  = math.sqrt(2*enerIIRight/mp)
        velDelII    = velIIRight - velIILeft
        for jj in loopTheta:
            thetaJJ = thetaTT[jj]
            for kk in loopPhi:
                phiKK    =  phiTT[kk]
                denTT    += deltaTheta*deltaPhi*velIISqr*velDelII*math.sin(thetaJJ)*distTT[ii,jj,kk]
                nByVxTT  += deltaTheta*deltaPhi*velIICube*velDelII*(-1)*(math.sin(thetaJJ))**2*math.cos(phiKK)*distTT[ii,jj,kk]
                nByVyTT  += deltaTheta*deltaPhi*velIICube*velDelII*(-1)*(math.sin(thetaJJ))**2*math.sin(phiKK)*distTT[ii,jj,kk]
                nByVzTT  += deltaTheta*deltaPhi*velIICube*velDelII*(-1)*math.sin(thetaJJ)*math.cos(thetaJJ)*distTT[ii,jj,kk]
                pxxItem1 += mp*deltaTheta*deltaPhi*velIIQuar*velDelII*(math.sin(thetaJJ))**3*(math.cos(phiKK))**2*distTT[ii,jj,kk]
                pyyItem1 += mp*deltaTheta*deltaPhi*velIIQuar*velDelII*(math.sin(thetaJJ))**3*(math.sin(phiKK))**2*distTT[ii,jj,kk]
                pzzItem1 += mp*deltaTheta*deltaPhi*velIIQuar*velDelII*math.sin(thetaJJ)*(math.cos(thetaJJ))**2*distTT[ii,jj,kk]


    den      = np.append(den, denTT)

    vxTT     = nByVxTT/denTT
    vyTT     = nByVyTT/denTT
    vzTT     = nByVzTT/denTT
    vtotalTT = math.sqrt(vxTT**2 + vyTT**2 + vzTT**2)
    vx       = np.append(vx, vxTT)
    vy       = np.append(vy, vyTT)
    vz       = np.append(vz, vzTT)
    vtotal   = np.append(vtotal, vtotalTT)

    pxxTT    = pxxItem1 - mp*denTT*vxTT**2
    pyyTT    = pyyItem1 - mp*denTT*vyTT**2
    pzzTT    = pzzItem1 - mp*denTT*vzTT**2
    pTT      = (pxxTT + pyyTT + pzzTT)/3
    pth      = np.append(pth, pTT)

densityM  = den###number density in m**(-3)
densityCM = den*10**(-6)###number density in cm*8(-3)
vxMS      = vx###velocity in m/s
vyMS      = vy###velocity in m/s
vzMS      = vz###velocity in m/s
vtotalMS  = vtotal###velocity in m/s
vxKMS     = vxMS*10**(-3)###velocity in km/s
vyKMS     = vyMS*10**(-3)###velocity in km/s
vzKMS     = vzMS*10**(-3)###velocity in km/s
vtotalKMS = vtotalMS*10**(-3)###velocity in km/s
pthPa     = pth###p in Pa
pthnPa    = pthPa*10**9###p in nPa

epochToPlot = epoch[timeArrToCal]





fileNameFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2016/12/17/mms1_fgm_brst_l2_20161217125504_v5.87.0.cdf')
infoFgm     = fileNameFgm.cdf_info()
epochFgm    = fileNameFgm.varget('Epoch')
bGSEVec     = fileNameFgm.varget('mms1_fgm_b_gse_brst_l2')

timeFgmTuple1    = np.where(epochFgm > epochStartType2000)
timeFgmTuple2    = np.where(epochFgm < epochStopType2000)
timeFgmArr1      = timeFgmTuple1[0]
timeFgmArr2      = timeFgmTuple2[0]
timeFgmArrToPlot = np.intersect1d(timeFgmArr1, timeFgmArr2)
epochFgmToPlot   = epochFgm[timeFgmArrToPlot]

btArr  = bGSEVec[:,3]
btArrToPlot  = btArr[timeFgmArrToPlot]

mu0  = 4*math.pi*10**(-7)
pb    = (btArrToPlot*10**(-9))**2/(2*mu0)
pbnPa = pb*10**9###magnetic pressure in nPa


fig, ax = plt.subplots()
xlim1 = [min(epochToPlot), max(epochToPlot)]
ax.plot(epochToPlot, pthnPa,color = 'darkorange')
ax.plot(epochFgmToPlot, pbnPa,color = 'deeppink')
ax.set_xlim(xlim1)
ax.set_ylabel('p [nPa]', fontsize = 18)
ax.set_ylim([0,0.6])


xtick1 = [2016,12,17,12,55, 10,000]
xtick2 = [2016,12,17,12,55, 20,000]
xtick3 = [2016,12,17,12,55, 30,000]
xtick4 = [2016,12,17,12,55, 40,000]

xtick1 = cdflib.cdfepoch.compute_tt2000(xtick1)
xtick2 = cdflib.cdfepoch.compute_tt2000(xtick2)
xtick3 = cdflib.cdfepoch.compute_tt2000(xtick3)
xtick4 = cdflib.cdfepoch.compute_tt2000(xtick4)

epochTick = [xtick1, xtick2, xtick3, xtick4]
strTick = ['12:55:10', '12:55:20', '12:55:30', '12:55:40']
ax.set_xticks(epochTick)
ax.set_xticklabels(strTick, fontsize = 15)


plt.show()




pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))

# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 09:48:28 2019
Modified on Mon Apr 8 16:50:23 2019

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
#import sympy
#import scipy.integrate
import matplotlib.pyplot as plt
###
import time


pyStartTime = time.time()
###
###constants==================
mp   = 1.67262158*10**(-27)
eV2J = 1.602176462*10**(-19)
kb   = 1.3806503*10**(-23)


timeInput1 = [2017,11,16,12,12,55,000]
timeInput2 = [2017,11,16,12,14,00,000]


fileName     = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2017/11/16/mms1_fpi_brst_l2_dis-dist_20171116121223_v3.3.0.cdf')
info         = fileName.cdf_info()
epoch        = fileName.varget('epoch')
epochPlusVar = fileName.varget('Epoch_plus_var')
#
startDelPhi = fileName.varget('mms1_dis_startdelphi_angle_brst')
#
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
dist0      = dist*10**12###dist in m s
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

    loopEner  = np.arange(0, enerTT.size)
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
                nByVxTT  += deltaTheta*deltaPhi*velIICube*velDelII*(math.sin(thetaJJ))**2*math.cos(phiKK)*distTT[ii,jj,kk]
                nByVyTT  += deltaTheta*deltaPhi*velIICube*velDelII*(math.sin(thetaJJ))**2*math.sin(phiKK)*distTT[ii,jj,kk]
                nByVzTT  += deltaTheta*deltaPhi*velIICube*velDelII*math.sin(thetaJJ)*math.cos(thetaJJ)*distTT[ii,jj,kk]
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
#fig, (ax1, ax2, ax3,ax4) = plt.subplots(4,1)


fig = plt.figure(figsize=(16,20),dpi=50)
fig.subplots_adjust(top = 0.90)
fig.subplots_adjust(bottom = 0.07)
fig.subplots_adjust(left = 0.13)
fig.subplots_adjust(right = 0.93)
ax1 = plt.subplot(8,1,1)
ax2 = plt.subplot(8,1,2)
ax3 = plt.subplot(8,1,3)
ax4 = plt.subplot(8,1,4)
ax5 = plt.subplot(8,1,5)
ax6 = plt.subplot(8,1,6)
ax7 = plt.subplot(8,1,7)
ax8 = plt.subplot(8,1,8)
#fig, (ax1, ax2, ax3,ax4) = plt.subplots(4,1)


xlim1 = [min(epochToPlot), max(epochToPlot)]
ax1.plot(epochToPlot, densityCM, 'g')
ax1.set_ylabel('Density wmm', fontsize = 10)
ax1.set_ylim([0,20])
ax1.set_xlim(xlim1)
ax1.set_title('Moments VS Official Moments', fontsize = 20)



ax3.plot(epochToPlot, vxKMS*(1), 'b')
ax3.plot(epochToPlot, vyKMS*(1), 'g')
ax3.plot(epochToPlot, vzKMS*(1), 'r')
ax3.set_ylabel('Velocity wmm', fontsize =10)
ax3.set_xlim(xlim1)
#ax3.set_ylim([-600, 200])



ax5.plot(epochToPlot, vtotalKMS, 'k')
ax5.set_ylabel('V total wmm', fontsize = 10)
ax5.set_xlim(xlim1)
ax5.set_ylim([0,600])


ax7.plot(epochToPlot, pthnPa, 'm')
ax7.set_ylabel('Thermal Pressure wmm', fontsize = 10)
ax7.set_xlim(xlim1)
ax7.set_ylim([0,0.4])




'''
part 2===============================================================================================================================
'''
momFileName = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-moms/2017/11/16/mms1_fpi_brst_l2_dis-moms_20171116121223_v3.3.0.cdf')
momInfo     = momFileName.cdf_info()
epochMom    = momFileName.varget('Epoch')
denMom      = momFileName.varget('mms1_dis_numberdensity_brst')
bulkvDBCS   = momFileName.varget('mms1_dis_bulkv_dbcs_brst')
tempPara    = momFileName.varget('mms1_dis_temppara_brst')
tempPerp    = momFileName.varget('mms1_dis_tempperp_brst')


timeMomTuple1         = np.where(epochMom > epochStartType2000)###This is tuple data
timeMomTuple2         = np.where(epochMom < epochStopType2000)###This is tuple data
timeMomArr1           = timeMomTuple1[0]
timeMomArr2           = timeMomTuple2[0]
timeMomArrToPlot       = np.intersect1d(timeMomArr1,timeMomArr2)

epochMomToPlot  = epochMom[timeMomArrToPlot]
denToPlot       = denMom[timeMomArrToPlot]
bulkvDBCSToPlot = bulkvDBCS[timeMomArrToPlot]
tempParaToPlot  = tempPara[timeMomArrToPlot]
tempPerpToPlot  = tempPerp[timeMomArrToPlot]

vxToPlot = bulkvDBCSToPlot[:,0]
vyToPlot = bulkvDBCSToPlot[:,1]
vzToPlot = bulkvDBCSToPlot[:,2]
vtotalToPlot = np.power((np.power(vxToPlot,2) + np.power(vyToPlot,2) + np.power(vzToPlot,2)), 1/2)

tempTotalToPlot  = 1/3*tempParaToPlot + 2/3*tempPerpToPlot
tempTotalToPlotJ = tempTotalToPlot * eV2J
denToPlotM       = denToPlot * 10**6
pthToPlotPa      = denToPlotM * tempTotalToPlotJ
pthToPlotnPa     = pthToPlotPa * 10**9




xlim2 = [min(epochMomToPlot), max(epochMomToPlot)]
ax2.plot(epochMomToPlot, denToPlot, 'g')
ax2.set_ylabel('Density official', fontsize = 10)
ax2.set_xlim(xlim2)
ax2.set_ylim([0,20])

ax4.plot(epochMomToPlot, vxToPlot, 'b')
ax4.plot(epochMomToPlot, vyToPlot, 'g')
ax4.plot(epochMomToPlot, vzToPlot, 'r')
ax4.set_ylabel('Velocity official', fontsize = 10)
ax4.set_xlim(xlim2)
ax4.set_ylim([-600,200])

ax6.plot(epochMomToPlot, vtotalToPlot, 'k')
ax6.set_ylabel('V total official', fontsize = 10)
ax6.set_xlim(xlim2)
ax6.set_ylim([0, 600])

ax8.plot(epochMomToPlot, pthToPlotnPa, 'm')
ax8.set_ylabel('Thermal Pressure official', fontsize = 10)
ax8.set_xlim(xlim2)
ax8.set_ylim([0,0.4])

plt.savefig('D:/shfaFigure/20171116121300/moments2.png')
plt.show()



pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))

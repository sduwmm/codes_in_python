# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 16:55:49 2019

@author: WORK
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 09:48:28 2019
Modified on Tue Apr 2 16:56:30 2019

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
import pyquaternion  ###better
from scipy import interpolate
###
import time

pyStartTime = time.time()
###
###constants==================
mp = 1.67262158 * 10 ** (-27)
eV2J = 1.602176462 * 10 ** (-19)
kb = 1.3806503 * 10 ** (-23)

timeInput1 = [2017, 11, 16, 12, 13, 00, 000]
timeInput2 = [2017, 11, 16, 12, 13, 10, 000]

epochStartType2000 = cdflib.cdfepoch.compute_tt2000(timeInput1)
epochStopType2000 = cdflib.cdfepoch.compute_tt2000(timeInput2)

'''
part 1===========================================================================================================================
'''

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

timeTuple1 = np.where(epoch > epochStartType2000)  ###This is tuple data
timeTuple2 = np.where(epoch < epochStopType2000)  ###This is tuple data
timeArr1 = timeTuple1[0]
timeArr2 = timeTuple2[0]
timeArrToCal = np.intersect1d(timeArr1, timeArr2)

###convert data to new data in SI Units
phi = phi * math.pi / 180  ###angle in radian
theta = theta * math.pi / 180  ###angle in radian
deltaPhi = phiDelta[0] * 2 * math.pi / 180  ###angle in radian
deltaTheta = thetaDelta[0] * 2 * math.pi / 180  ###angle in radian
dist0 = dist * 10 ** 12  ###dist in m s
energyJ = energy * eV2J  ###energy in J
enerDeltaJ = energyDelta * eV2J  ###energy in J
###convert data to new data in SI Units


den = []
vx = []
vy = []
vz = []
vtotal = []
pth = []
for tt in timeArrToCal:
    phiTT = phi[tt]
    thetaTT = theta
    distTT = dist0[tt]
    enerTT = energyJ[tt]
    enerDeltaTT = enerDeltaJ[tt]

    loopEner = np.arange(0, enerTT.size)
    loopPhi = np.arange(0, phiTT.size)
    loopTheta = np.arange(0, thetaTT.size)

    denTT = 0
    nByVxTT = 0
    nByVyTT = 0
    nByVzTT = 0
    pxxItem1 = 0
    pyyItem1 = 0
    pzzItem1 = 0
    for ii in loopEner:
        enerII = enerTT[ii]
        enerDeltaII = enerDeltaTT[ii]
        velII = math.sqrt(2 * enerII / mp)
        velIISqr = velII ** 2
        velIICube = velII ** 3
        velIIQuar = velII ** 4
        enerIILeft = enerII - enerDeltaII
        enerIIRight = enerII + enerDeltaII
        velIILeft = math.sqrt(2 * enerIILeft / mp)
        velIIRight = math.sqrt(2 * enerIIRight / mp)
        velDelII = velIIRight - velIILeft
        for jj in loopTheta:
            thetaJJ = thetaTT[jj]
            for kk in loopPhi:
                phiKK = phiTT[kk]
                denTT += deltaTheta * deltaPhi * velIISqr * velDelII * math.sin(thetaJJ) * distTT[ii, jj, kk]
                nByVxTT += (-1)*deltaTheta * deltaPhi * velIICube * velDelII * (math.sin(thetaJJ)) ** 2 * math.cos(phiKK) * \
                           distTT[ii, jj, kk]
                nByVyTT += (-1)*deltaTheta * deltaPhi * velIICube * velDelII * (math.sin(thetaJJ)) ** 2 * math.sin(phiKK) * \
                           distTT[ii, jj, kk]
                nByVzTT += (-1)*deltaTheta * deltaPhi * velIICube * velDelII * math.sin(thetaJJ) * math.cos(thetaJJ) * \
                           distTT[ii, jj, kk]
                pxxItem1 += mp * deltaTheta * deltaPhi * velIIQuar * velDelII * (math.sin(thetaJJ)) ** 3 * (
                    math.cos(phiKK)) ** 2 * distTT[ii, jj, kk]
                pyyItem1 += mp * deltaTheta * deltaPhi * velIIQuar * velDelII * (math.sin(thetaJJ)) ** 3 * (
                    math.sin(phiKK)) ** 2 * distTT[ii, jj, kk]
                pzzItem1 += mp * deltaTheta * deltaPhi * velIIQuar * velDelII * math.sin(thetaJJ) * (
                    math.cos(thetaJJ)) ** 2 * distTT[ii, jj, kk]

    den = np.append(den, denTT)

    vxTT = nByVxTT / denTT
    vyTT = nByVyTT / denTT
    vzTT = nByVzTT / denTT
    vtotalTT = math.sqrt(vxTT ** 2 + vyTT ** 2 + vzTT ** 2)
    vx = np.append(vx, vxTT)
    vy = np.append(vy, vyTT)
    vz = np.append(vz, vzTT)
    vtotal = np.append(vtotal, vtotalTT)

    pxxTT = pxxItem1 - mp * denTT * vxTT ** 2
    pyyTT = pyyItem1 - mp * denTT * vyTT ** 2
    pzzTT = pzzItem1 - mp * denTT * vzTT ** 2
    pTT = (pxxTT + pyyTT + pzzTT) / 3
    pth = np.append(pth, pTT)

densityM = den  ###number density in m**(-3)
densityCM = den * 10 ** (-6)  ###number density in cm*8(-3)
vxMS = vx  ###velocity in m/s
vyMS = vy  ###velocity in m/s
vzMS = vz  ###velocity in m/s
vtotalMS = vtotal  ###velocity in m/s
vxKMS = vxMS * 10 ** (-3)  ###velocity in km/s
vyKMS = vyMS * 10 ** (-3)  ###velocity in km/s
vzKMS = vzMS * 10 ** (-3)  ###velocity in km/s
vtotalKMS = vtotalMS * 10 ** (-3)  ###velocity in km/s
pthPa = pth  ###p in Pa
pthnPa = pthPa * 10 ** 9  ###p in nPa

epochToPlot = epoch[timeArrToCal]

# fig, (ax1, ax2, ax3,ax4) = plt.subplots(4,1)

'''
part 2===============================================================================================================================
'''
mecFileName = cdflib.CDF(
    'D:/data/mms/mms1/mec/brst/l2/ephts04d/2017/11/16/mms1_mec_brst_l2_ephts04d_20171116121225_v2.2.0.cdf')
mecInfo = mecFileName.cdf_info()
mecEpoch = mecFileName.varget('Epoch')
quatECI2DBCSArr = mecFileName.varget('mms1_mec_quat_eci_to_dbcs')
quatECI2GSEArr = mecFileName.varget('mms1_mec_quat_eci_to_gse')

timeMecTuple1 = np.where(mecEpoch > epochStartType2000)  ###This is tuple data
timeMecTuple2 = np.where(mecEpoch < epochStopType2000)  ###This is tuple data
timeMecArr1 = timeMecTuple1[0]
timeMecArr2 = timeMecTuple2[0]
timeMecArrToCal = np.intersect1d(timeMecArr1, timeMecArr2)

# timeMecArrToCal = np.delete(timeMecArrToCal, 2165, axis = 0)
# timeMecArrToCal = np.delete(timeMecArrToCal, 0,    axis = 0)


quatECI2DBCSArrToCal = quatECI2DBCSArr[timeMecArrToCal]
quatECI2GSEArrToCal  = quatECI2GSEArr[timeMecArrToCal]
mecEpochToPlot       = mecEpoch[timeMecArrToCal]
# mecEpochToPlot[2165] = max(epochToPlot)
mecEpochToPlot[0] = min(epochToPlot)

f1 = interpolate.interp1d(mecEpochToPlot, quatECI2DBCSArrToCal[:, 0])
quatECI2DBCSArrToCalW = f1(epochToPlot)
quatECI2DBCSArrToCalW = quatECI2DBCSArrToCalW.reshape(-1,1)
f1 = interpolate.interp1d(mecEpochToPlot, quatECI2DBCSArrToCal[:, 1])
quatECI2DBCSArrToCalX = f1(epochToPlot)
quatECI2DBCSArrToCalX = quatECI2DBCSArrToCalX.reshape(-1,1)
f1 = interpolate.interp1d(mecEpochToPlot, quatECI2DBCSArrToCal[:, 2])
quatECI2DBCSArrToCalY = f1(epochToPlot)
quatECI2DBCSArrToCalY = quatECI2DBCSArrToCalY.reshape(-1,1)
f1 = interpolate.interp1d(mecEpochToPlot, quatECI2DBCSArrToCal[:, 3])
quatECI2DBCSArrToCalZ = f1(epochToPlot)
quatECI2DBCSArrToCalZ = quatECI2DBCSArrToCalZ.reshape(-1,1)
quatECI2DBCSArrToCalNew = np.hstack((quatECI2DBCSArrToCalW, quatECI2DBCSArrToCalX, quatECI2DBCSArrToCalY, quatECI2DBCSArrToCalZ))


f2 = interpolate.interp1d(mecEpochToPlot, quatECI2GSEArrToCal[:, 0])
quatECI2GSEArrToCalW = f2(epochToPlot)
quatECI2GSEArrToCalW = quatECI2GSEArrToCalW.reshape(-1,1)
f2 = interpolate.interp1d(mecEpochToPlot, quatECI2GSEArrToCal[:, 1])
quatEIC2GSEArrToCalX = f2(epochToPlot)
quatEIC2GSEArrToCalX = quatEIC2GSEArrToCalX.reshape(-1,1)
f2 = interpolate.interp1d(mecEpochToPlot, quatECI2GSEArrToCal[:, 2])
quatEIC2GSEArrToCalY = f2(epochToPlot)
quatEIC2GSEArrToCalY = quatEIC2GSEArrToCalY.reshape(-1,1)
f2 = interpolate.interp1d(mecEpochToPlot, quatECI2GSEArrToCal[:, 3])
quatEIC2GSEArrToCalZ = f2(epochToPlot)
quatEIC2GSEArrToCalZ = quatEIC2GSEArrToCalZ.reshape(-1,1)
quatECI2GSEArrToCalNew = np.hstack((quatECI2GSEArrToCalW, quatEIC2GSEArrToCalX, quatEIC2GSEArrToCalY, quatEIC2GSEArrToCalZ))


vGSE = [0, 0, 0]

for tt in np.arange(quatECI2DBCSArrToCalW.size):
    print(tt)
    vxDBCSTT = vxKMS[tt]
    vyDBCSTT = vyKMS[tt]
    vzDBCSTT = vzKMS[tt]
    vDBCSTT = np.hstack((vxDBCSTT, vyDBCSTT, vzDBCSTT))
    quatECI2DBCSTT = pyquaternion.Quaternion(quatECI2DBCSArrToCalNew[tt, :])
    quatECI2DBCSTT = quatECI2DBCSTT.normalised
    quatDBCS2ECITT = quatECI2DBCSTT.inverse
    vECITT = quatDBCS2ECITT.rotate(vDBCSTT)
    quatECI2GSETT = pyquaternion.Quaternion(quatECI2GSEArrToCalNew[tt, :])
    quatECI2GSETT = quatECI2GSETT.normalised
    vGSETT = quatECI2GSETT.rotate(vECITT)
    vGSE = np.vstack((vGSE, vGSETT))  ###???

vGSE = np.delete(vGSE, 0, axis=0)

'''
================================================================================================================================part2
'''

'''
part 3===============================================================================================================================
'''
momFileName = cdflib.CDF(
    'D:/data/mms/mms1/fpi/brst/l2/dis-moms/2017/11/16/mms1_fpi_brst_l2_dis-moms_20171116121223_v3.3.0.cdf')
momInfo = momFileName.cdf_info()
epochMom = momFileName.varget('Epoch')
denMom = momFileName.varget('mms1_dis_numberdensity_brst')
bulkvDBCS = momFileName.varget('mms1_dis_bulkv_dbcs_brst')
bulkvGSE  = momFileName.varget('mms1_dis_bulkv_gse_brst')
tempPara = momFileName.varget('mms1_dis_temppara_brst')
tempPerp = momFileName.varget('mms1_dis_tempperp_brst')

timeMomTuple1 = np.where(epochMom > epochStartType2000)  ###This is tuple data
timeMomTuple2 = np.where(epochMom < epochStopType2000)  ###This is tuple data
timeMomArr1 = timeMomTuple1[0]
timeMomArr2 = timeMomTuple2[0]
timeMomArrToPlot = np.intersect1d(timeMomArr1, timeMomArr2)

epochMomToPlot = epochMom[timeMomArrToPlot]
denToPlot = denMom[timeMomArrToPlot]
bulkvDBCSToPlot = bulkvDBCS[timeMomArrToPlot]
bulkvGSEToPlot  = bulkvGSE[timeMomArrToPlot]
tempParaToPlot = tempPara[timeMomArrToPlot]
tempPerpToPlot = tempPerp[timeMomArrToPlot]

vxToPlot = bulkvDBCSToPlot[:, 0]
vyToPlot = bulkvDBCSToPlot[:, 1]
vzToPlot = bulkvDBCSToPlot[:, 2]
vtotalToPlot = np.power((np.power(vxToPlot, 2) + np.power(vyToPlot, 2) + np.power(vzToPlot, 2)), 1 / 2)

tempTotalToPlot = 1 / 3 * tempParaToPlot + 2 / 3 * tempPerpToPlot
tempTotalToPlotJ = tempTotalToPlot * eV2J
denToPlotM = denToPlot * 10 ** 6
pthToPlotPa = denToPlotM * tempTotalToPlotJ
pthToPlotnPa = pthToPlotPa * 10 ** 9

fig = plt.figure(figsize=(16, 20), dpi=50)
fig.subplots_adjust(top=0.90)
fig.subplots_adjust(bottom=0.07)
fig.subplots_adjust(left=0.13)
fig.subplots_adjust(right=0.93)
ax1 = plt.subplot(6, 1, 1)
ax2 = plt.subplot(6, 1, 2)
ax3 = plt.subplot(6, 1, 3)
ax4 = plt.subplot(6, 1, 4)
ax5 = plt.subplot(6, 1, 5)
ax6 = plt.subplot(6, 1, 6)

# fig, (ax1, ax2, ax3,ax4) = plt.subplots(4,1)


xlim1 = [min(epochToPlot), max(epochToPlot)]
ax1.plot(epochToPlot, vxKMS, 'darkorange')
ax1.plot(epochToPlot, bulkvDBCSToPlot[:,0], 'deepskyblue')
ax1.set_ylabel('Vx DBCS', fontsize=10)
ax1.set_xlim(xlim1)
ax1.set_ylim([-500, 0])
ax1.set_title('Moments VS Official Moments', fontsize=20)


ax2.plot(epochToPlot, vyKMS, 'darkorange')
ax2.plot(epochToPlot, bulkvDBCSToPlot[:,1], 'deepskyblue')
ax2.set_xlim(xlim1)
ax2.set_ylim([0,500])
ax2.set_ylabel('Vy DBCS', fontsize=10)

ax3.plot(epochToPlot, vzKMS, 'darkorange')
ax3.plot(epochToPlot, bulkvDBCSToPlot[:,2], 'deepskyblue')
ax3.set_ylabel('Vz DBCS', fontsize = 10)
ax3.set_xlim(xlim1)
ax3.set_ylim([-250, 250])

ax4.plot(epochToPlot, vGSE[:,0], 'darkorange')
ax4.plot(epochToPlot, bulkvGSEToPlot[:,0], 'deepskyblue')
ax4.set_xlim(xlim1)
ax4.set_ylabel('Vx GSE', fontsize = 10)
ax4.set_ylim([-500, 0])

ax5.plot(epochToPlot, vGSE[:,1], 'darkorange')
ax5.plot(epochToPlot, bulkvGSEToPlot[:,1], 'deepskyblue')
ax5.set_xlim(xlim1)
ax5.set_ylabel('Vy GSE', fontsize = 10)
ax5.set_ylim([0,500])

ax6.plot(epochToPlot, vGSE[:,2], 'darkorange')
ax6.plot(epochToPlot, bulkvGSEToPlot[:,2], 'deepskyblue')
ax6.set_xlim(xlim1)
ax6.set_ylabel('Vz GSE', fontsize = 10)
ax6.set_ylim([-250, 250])


plt.savefig('D:/shfaFigure/20171116121300/moments/dbcs_gse.png')
plt.show()

pyEndTime = time.time()
print(
    '==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))

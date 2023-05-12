'''
Created on Tue Apr 9 11:09:11 2019
Changed on Mon Mar 22 14:19:00 2021
@author: Mengmeng Wang at sdu
'''

'''
calculate the zero-, first, second- moments 
of MMS FPI ion distribution data
from the phase space density
in given time,
phi angle and theta angle range,
energy range

solar wind beam

Version 0 



!!!Note
dist: time*energy*theta*phi=epoch_number*32*16*32
'''
'''
The certain energy range, theta range, phi range are inputted 
for moments excluding solar wid ions.
'''

import cdflib
import math
import numpy as np
import matplotlib.pyplot as plt
import time

pyStartTime = time.time()
def MomentsNoSWDuration(Duration, FileDist, Theta, Phi):
    ###constants==================
    mp = 1.67262158 * 10 ** (-27)
    eV2J = 1.602176462 * 10 ** (-19)
    kb = 1.3806503 * 10 ** (-23)
    mu0 = 4 * math.pi * 10 ** (-7)
    ###==================constants
    ###init====================
    #EnerSW = np.array(Ener)
    PhiSW = np.array(Phi)
    ThetaSW = np.array(Theta)
    ###========================init
    ###get data==========================
    InfoDist = FileDist.cdf_info()
    Epoch = FileDist.varget('Epoch')
    DistRaw = FileDist.varget('mms2_dis_dist_brst')
    DistErrRaw = FileDist.varget('mms2_dis_disterr_brst')
    EnergyRaw = FileDist.varget('mms2_dis_energy_brst')
    EnergyDeltaRaw = FileDist.varget('mms2_dis_energy_delta_brst')
    ThetaRaw = FileDist.varget('mms2_dis_theta_brst')
    ThetaDeltaRaw = FileDist.varget('mms2_dis_theta_delta_brst')
    PhiRaw = FileDist.varget('mms2_dis_phi_brst')
    PhiDeltaRaw = FileDist.varget('mms2_dis_phi_delta_brst')
    ###==========================get data
    ###Get the epoch slice====================================
    EpochStartType2000 = cdflib.cdfepoch.compute_tt2000(Duration[0])
    EpochEndType2000 = cdflib.cdfepoch.compute_tt2000(Duration[1])
    TimeTplStart = np.where(Epoch >= EpochStartType2000)  ###This is tuple data
    TimeTplEnd = np.where(Epoch <= EpochEndType2000)  ###This is tuple data
    TimeArr1 = TimeTplStart[0]
    TimeArr2 = TimeTplEnd[0]
    TimeArrToCal = np.intersect1d(TimeArr1, TimeArr2)
    ####====================Get the epoch slice
    ###convert data to new data in SI Units
    dist0 = DistRaw * 10 ** 12  ###dist in s3/m6
    distErr = DistErrRaw * 10 ** 12  ### dist in s3/m6
    energyJ = EnergyRaw * eV2J  ###energy in J
    enerDeltaJ = EnergyDeltaRaw * eV2J  ###energy in J
    deltaPhi = PhiDeltaRaw[0] * 2 * math.pi / 180  ###angle in radian
    deltaTheta = ThetaDeltaRaw[0] * 2 * math.pi / 180  ###angle in radian
    phi = PhiRaw * math.pi / 180  ###angle in radian
    theta = ThetaRaw * math.pi / 180  ###angle in radian
    ###convert data to new data in SI Units

    den = []
    vx = []
    vy = []
    vz = []
    vtotal = []
    pth = []
    for tt in TimeArrToCal:
        phiTT = phi[tt]
        thetaTT = theta
        distTT = dist0[tt]
        distErrTT = distErr[tt]###inserted by wmm on 2021-03-22
        enerTT = energyJ[tt]
        enerDeltaTT = enerDeltaJ[tt]
        distTT[distTT < distErrTT] = 0###inserted by wmm on 2021-03-22

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
                    if jj in ThetaSW:  ###added by wmm 2021-03-16
                        if kk in PhiSW:  ###added by wmm 2021-03-16
                            distTT[:, jj, kk] = 0  ###added by wmm 2021-03-16
                    denTT += deltaTheta * deltaPhi * velIISqr * velDelII * math.sin(thetaJJ) * distTT[ii, jj, kk]
                    nByVxTT += deltaTheta * deltaPhi * velIICube * velDelII * (-1) * (
                        math.sin(thetaJJ)) ** 2 * math.cos(phiKK) * distTT[ii, jj, kk]
                    nByVyTT += deltaTheta * deltaPhi * velIICube * velDelII * (-1) * (
                        math.sin(thetaJJ)) ** 2 * math.sin(phiKK) * distTT[ii, jj, kk]
                    nByVzTT += deltaTheta * deltaPhi * velIICube * velDelII * (-1) * math.sin(thetaJJ) * math.cos(
                        thetaJJ) * distTT[ii, jj, kk]
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

    return (densityCM, vxKMS, vyKMS, vzKMS, vtotalKMS, pthnPa, TimeArrToCal)

def MomentsSWDuration(Duration, FileDist, Energy, Theta, Phi):
    ###constants==================
    mp = 1.67262158 * 10 ** (-27)
    eV2J = 1.602176462 * 10 ** (-19)
    kb = 1.3806503 * 10 ** (-23)
    mu0 = 4 * math.pi * 10 ** (-7)
    ###==================constants
    ###init====================
    #EnerSW = np.array(Ener)
    PhiSW = np.array(Phi)
    ThetaSW = np.array(Theta)
    ###========================init
    ###get data==========================
    InfoDist = FileDist.cdf_info()
    Epoch = FileDist.varget('Epoch')
    DistRaw = FileDist.varget('mms2_dis_dist_brst')
    DistErrRaw = FileDist.varget('mms2_dis_disterr_brst')
    EnergyRaw = FileDist.varget('mms2_dis_energy_brst')
    EnergyDeltaRaw = FileDist.varget('mms2_dis_energy_delta_brst')
    ThetaRaw = FileDist.varget('mms2_dis_theta_brst')
    ThetaDeltaRaw = FileDist.varget('mms2_dis_theta_delta_brst')
    PhiRaw = FileDist.varget('mms2_dis_phi_brst')
    PhiDeltaRaw = FileDist.varget('mms2_dis_phi_delta_brst')
    ###==========================get data
    ###Get the epoch slice====================================
    EpochStartType2000 = cdflib.cdfepoch.compute_tt2000(Duration[0])
    EpochEndType2000 = cdflib.cdfepoch.compute_tt2000(Duration[1])
    TimeTplStart = np.where(Epoch >= EpochStartType2000)  ###This is tuple data
    TimeTplEnd = np.where(Epoch <= EpochEndType2000)  ###This is tuple data
    TimeArr1 = TimeTplStart[0]
    TimeArr2 = TimeTplEnd[0]
    TimeArrToCal = np.intersect1d(TimeArr1, TimeArr2)
    ####====================Get the epoch slice
    ###convert data to new data in SI Units
    dist0 = DistRaw * 10 ** 12  ###dist in s3/m6
    distErr = DistErrRaw * 10 ** 12  ### dist in s3/m6
    energyJ = EnergyRaw * eV2J  ###energy in J
    enerDeltaJ = EnergyDeltaRaw * eV2J  ###energy in J
    deltaPhi = PhiDeltaRaw[0] * 2 * math.pi / 180  ###angle in radian
    deltaTheta = ThetaDeltaRaw[0] * 2 * math.pi / 180  ###angle in radian
    phi = PhiRaw * math.pi / 180  ###angle in radian
    theta = ThetaRaw * math.pi / 180  ###angle in radian
    ###convert data to new data in SI Units

    den = []
    vx = []
    vy = []
    vz = []
    vtotal = []
    pth = []
    for tt in TimeArrToCal:
        phiTT = phi[tt]
        thetaTT = theta
        distTT = dist0[tt]
        distErrTT = distErr[tt]###inserted by wmm on 2021-03-22
        enerTT = energyJ[tt]
        enerDeltaTT = enerDeltaJ[tt]
        distTT[distTT < distErrTT] = 0###inserted by wmm on 2021-03-22

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
                    if jj in ThetaSW:  ###added by wmm 2021-03-16
                        if kk in PhiSW:  ###added by wmm 2021-03-16
                            denTT += deltaTheta * deltaPhi * velIISqr * velDelII * math.sin(thetaJJ) * distTT[ii, jj, kk]
                            nByVxTT += deltaTheta * deltaPhi * velIICube * velDelII * (-1) * (
                                math.sin(thetaJJ)) ** 2 * math.cos(phiKK) * distTT[ii, jj, kk]
                            nByVyTT += deltaTheta * deltaPhi * velIICube * velDelII * (-1) * (
                                math.sin(thetaJJ)) ** 2 * math.sin(phiKK) * distTT[ii, jj, kk]
                            nByVzTT += deltaTheta * deltaPhi * velIICube * velDelII * (-1) * math.sin(thetaJJ) * math.cos(
                                thetaJJ) * distTT[ii, jj, kk]
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

    return (densityCM, vxKMS, vyKMS, vzKMS, vtotalKMS, pthnPa, TimeArrToCal)

TimeInputStart = [2019, 2, 23, 5, 22, 5, 921]
TimeInputEnd = [2019, 2, 23, 5, 22, 22, 0]
FileDistInput = \
    cdflib.CDF('D:/data/mms/mms2/fpi/brst/l2/dis-dist/2019/02/23/mms2_fpi_brst_l2_dis-dist_20190223052123_v3.3.0.cdf')
DurationInput = [TimeInputStart, TimeInputEnd]
EnergySWInput = np.arange(13, 21, 1)
ThetaSWInput = np.arange(3, 11, 1)
PhiSWInput = [0,1,2,3,29,30,31]
momentsFore = MomentsNoSWDuration(DurationInput, FileDistInput, ThetaSWInput, PhiSWInput)
momentsSW = MomentsSWDuration(DurationInput, FileDistInput, EnergySWInput, ThetaSWInput, PhiSWInput)
densityFore = momentsFore[0]
vxFore = momentsFore[1]
vyFore = momentsFore[2]
vzFore = momentsFore[3]
vtotalFore = momentsFore[4]
pthnPaFore = momentsFore[5]
TimeArrToCal = momentsFore[6]
Epoch = FileDistInput.varget('Epoch')
EpochPlot = Epoch[TimeArrToCal]

densitySW = momentsSW[0]
vxSW = momentsSW[1]
vySW = momentsSW[2]
vzSW = momentsSW[3]
vtotalSW = momentsSW[4]
pthnPaSW = momentsSW[5]
###=================================================================================
fig = plt.figure(figsize=(8,10),dpi=300)
fig.subplots_adjust(top = 0.90)
fig.subplots_adjust(bottom = 0.07)
fig.subplots_adjust(left = 0.13)
fig.subplots_adjust(right = 0.93)
ax1 = plt.subplot(6,1,1)
ax2 = plt.subplot(6,1,2)
ax3 = plt.subplot(6,1,3)
ax4 = plt.subplot(6,1,4)
ax5 = plt.subplot(6,1,5)
ax6 = plt.subplot(6,1,6)
##
xlim1 = [min(EpochPlot), max(EpochPlot)]
ax1.plot(EpochPlot, densityFore, 'deepskyblue')
ax1.plot(EpochPlot, densitySW, 'darkorange')
ax1.set_ylabel('Density', fontsize = 15)
ax1.set_ylim([0,15])
ax1.set_xlim(xlim1)
ax1.set_title('Foreshock transient Liu MMS2', fontsize = 20)
#
ax2.plot(EpochPlot, vxFore, 'deepskyblue')
ax2.plot(EpochPlot, vxSW, 'darkorange')
ax2.set_ylabel('Vx', fontsize =15)
ax2.set_xlim(xlim1)
ax3.set_ylim([-400, 100])
#
ax3.plot(EpochPlot, vyFore, 'deepskyblue')
ax3.plot(EpochPlot, vySW, 'darkorange')
ax3.set_ylabel('Vy', fontsize = 15)
ax3.set_xlim(xlim1)
ax3.set_ylim([-100,200])
#
ax4.plot(EpochPlot, vzFore, 'deepskyblue')
ax4.plot(EpochPlot, vzSW, 'darkorange')
ax4.set_ylabel('Vz', fontsize = 15)
ax4.set_xlim(xlim1)
ax4.set_ylim([-150,250])
#
ax5.plot(EpochPlot, vtotalFore, 'deepskyblue')
ax5.plot(EpochPlot, vtotalSW, 'darkorange')
ax5.set_ylabel('Vtotal', fontsize = 15)
ax5.set_xlim(xlim1)
ax5.set_ylim([0, 400])
#
ax6.plot(EpochPlot, pthnPaFore, 'deepskyblue')
ax6.plot(EpochPlot, pthnPaSW, 'darkorange')
ax6.set_ylabel('Thermal Pressure', fontsize = 15)
ax6.set_xlim(xlim1)
ax6.set_ylim([0,0.3])
#



plt.savefig('D:\\OneDrive - mail.sdu.edu.cn\\work\code&figure\\figure2021/MMS2_2019_.png')
plt.show()
plt.close()



pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))




'''
Created on Tue Apr 9 11:09:20 2019
Changed on Tue Mar 23 17:33:39 2021
Test on June 17 19:00:00 2021
Modified on June 20 17:15:20 2021
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
!!!Note
dist slice rules changed after I update the Python, numpy and PyCharm
Written on 2021-06-17
'''
'''
!!!Note
dist: time*phi*theta*energy=epoch_number*32*16*32
Written on 2021-06-17
['Epoch', 'Epoch_plus_var', 'Epoch_minus_var', 
'mms1_des_errorflags_brst', 'mms1_des_compressionloss_brst', 
'mms1_des_steptable_parity_brst', 'mms1_des_startdelphi_count_brst', 
'mms1_des_startdelphi_angle_brst', 
'mms1_des_phi_brst', 'mms1_des_phi_delta_brst', 
'mms1_des_dist_brst', 'mms1_des_disterr_brst', 
'mms1_des_sector_despinp_brst', 'mms1_des_sector_label_brst', 
'mms1_des_pixel_label_brst', 'mms1_des_energy_label_brst', 
'mms1_des_theta_brst', 'mms1_des_theta_delta_brst', 
'mms1_des_energy_brst', 'mms1_des_energy_delta_brst']
'''

import cdflib
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import time

###constants===============================
mp = 1.67262158 * 10 ** (-27)
eV2J = 1.602176462 * 10 ** (-19)
kb = 1.3806503 * 10 ** (-23)
mu0 = 4 * math.pi * 10 ** (-7)
###=================================constants

pyStartTime = time.time()
def SkyMap(Duration, FileDist, Theta, Phi, FigurePath):
    ###init============================
    ThetaNo=np.array(Theta)
    PhiNo = np.array(Phi)
    ###===================================init
    ###get data==========================
    InfoDist=FileDist.cdf_info()
    EpochData=FileDist.varget('Epoch')
    PhiData=FileDist.varget('mms1_des_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_des_phi_delta_brst')
    DistData=FileDist.varget('mms1_des_dist_brst')
    DistErrData=FileDist.varget('mms1_des_disterr_brst')
    ThetaData=FileDist.varget('mms1_des_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_des_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_des_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_des_energy_delta_brst')
    ###===================================get data
    ###Get the epoch slice==================================================
    EpochStartType2000=cdflib.cdfepoch.compute_tt2000(Duration[0])
    EpochEndType2000=cdflib.cdfepoch.compute_tt2000(Duration[1])
    EpochTplStart=np.where(EpochData >= EpochStartType2000)
    EpochTplEnd=np.where(EpochData <= EpochEndType2000)
    EpochArr1=EpochTplStart[0]
    EpochArr2=EpochTplEnd[0]
    EpochArrToCal=np.intersect1d(EpochArr1, EpochArr2)
    ####=================================Get the epoch slice

    for tt in EpochArrToCal:
        ###
        DistTT = DistData[tt]
        DistErrTT = DistErrData[tt]
        EnergyChannelTT = EnergyChannel[tt]
        EnergyChannelDeltaTT = EnergyChannelDelta[tt]
        # ThetaTT=ThetaData
        # PhiTT=PhiData[EpochInd]

        DistAlongEnergySum = np.zeros([32, 16])  ##changed on 2021-06-17

        #for ee in np.arange(11, 32, 1):  # without bins above 100V 2021-06-20
        #for ee in np.arange(4,5,1):#without bins 2021-06-20
        #for ee in np.arange(2, 32, 1):  # without bins below 5eV 2021-06-17
        for ee in np.arange(2, 32, 1):

            EnergyInd = ee
            DistTTEE = DistTT[:, :, ee]  # take out the dist of certain energy channel
            DistErrTTEE = DistErrTT[:, :, ee]  # take out the dist error of certain energy channel
            #

            for jj in PhiNo:
                for kk in ThetaNo:
                    if jj in PhiNo:
                        if kk in ThetaNo:
                            DistTTEE[jj, kk] = 0

            #
            #
            DistTTEE[DistTTEE <= DistErrTTEE] = 0
            #
            DistAlongEnergySum = DistAlongEnergySum + DistTTEE
            #
            DistSumLg = np.log10(DistAlongEnergySum)
            DistSumLg = DistSumLg.transpose()  ##added on 2021-06-17

        #
        Figure = plt.figure(figsize=(5, 4), dpi=100)
        AX = plt.subplot(1, 1, 1)
        CMap = AX.pcolor(DistSumLg, cmap='jet', vmin=-25.2, vmax=-23.9)
        CBar = Figure.colorbar(CMap, ax=AX)
        CBar.set_label(r'$Log_{10} PSD (s^3/cm^6)$', x=0, y=0.55, rotation=90, fontsize=10)
        #
        Time = Duration[0]
        TimeStr = str(Time[0]) + '-' + str(Time[1]) + '-' + str(Time[2]) + '/' + \
                  str(Time[3]) + ':' + str(Time[4]) + ':' + str(Time[5]) + '.' + \
                  str(Time[6]) + '+' + str(tt-EpochArrToCal[0])
        AX.text(4, 16.5, TimeStr, fontsize=10)
        if not os.path.exists(FigurePath):
            os.makedirs(FigurePath)
        FigureName = FigurePath + 'MMS1_' + str(Time[0]) + str(Time[1]) + str(Time[2]) + \
                     '-' + str(Time[3]) + str(Time[4]) + str(Time[5]) + '.' + str(Time[6]) + \
                    '_' + '+' + str(tt-EpochArrToCal[0]) + \
                     '.png'

        plt.savefig(FigureName)
        plt.show()
        plt.close()
    return()
    #return(DistTT20)


def SkyMapTrans(Duration, FileDist, Theta, Phi, FigurePath):
    ###init============================
    ThetaNo=np.array(Theta)
    PhiNo = np.array(Phi)
    ###===================================init
    ###get data==========================
    InfoDist=FileDist.cdf_info()
    EpochData=FileDist.varget('Epoch')
    PhiData=FileDist.varget('mms1_des_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_des_phi_delta_brst')
    DistData=FileDist.varget('mms1_des_dist_brst')
    DistErrData=FileDist.varget('mms1_des_disterr_brst')
    ThetaData=FileDist.varget('mms1_des_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_des_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_des_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_des_energy_delta_brst')
    ###===================================get data
    ###Get the epoch slice==================================================
    EpochStartType2000=cdflib.cdfepoch.compute_tt2000(Duration[0])
    EpochEndType2000=cdflib.cdfepoch.compute_tt2000(Duration[1])
    EpochTplStart=np.where(EpochData >= EpochStartType2000)
    EpochTplEnd=np.where(EpochData <= EpochEndType2000)
    EpochArr1=EpochTplStart[0]
    EpochArr2=EpochTplEnd[0]
    EpochArrToCal=np.intersect1d(EpochArr1, EpochArr2)
    ####=================================Get the epoch slice

    for tt in EpochArrToCal:####Epoch index
        DistTT = DistData[tt]
        DistErrTT = DistErrData[tt]
        EnergyChannelTT = EnergyChannel[tt]
        EnergyChannelDeltaTT = EnergyChannelDelta[tt]
        # ThetaTT=ThetaData
        # PhiTT=PhiData[EpochInd]

        DistAlongEnergySum = np.zeros([32, 16])  ##changed on 2021-06-17

        #for ee in np.arange(2, 32, 1):  # without bins below 5eV 2021-06-17
        for ee in np.arange(0, 32, 1):  

            EnergyInd = ee
            DistTTEE = DistTT[:, :, ee]  # take out the dist of certain energy channel
            DistErrTTEE = DistErrTT[:, :, ee]  # take out the dist error of certain energy channel
            #

            for jj in PhiNo:
                for kk in ThetaNo:
                    if jj in PhiNo:
                        if kk in ThetaNo:
                            DistTTEE[jj, kk] = 0

            #
            #
            DistTTEE[DistTTEE <= DistErrTTEE] = 0
            #
            DistAlongEnergySum = DistAlongEnergySum + DistTTEE
            #
            DistSumLg = np.log10(DistAlongEnergySum)
            DistSumLg = DistSumLg.transpose()  ##added on 2021-06-17
            #
            [DistA, DistB] = np.hsplit(DistSumLg, 2)
            DistFirst = np.hstack((DistB, DistA))
            DistSecond = np.zeros([16, 32])
            # reverse for list
            for dd in np.arange(0, 16, 1):
                DistSecond[dd] = DistFirst[15 - dd]

        #
        Figure = plt.figure(figsize=(5, 4), dpi=100)
        AX = plt.subplot(1, 1, 1)
        # CMap=AX.pcolor(DistSumLg,cmap='jet',vmin=-26,vmax=-20)
        CMap = AX.pcolor(DistSecond, cmap='jet', vmin=-26, vmax=-20)
        CBar = Figure.colorbar(CMap, ax=AX)
        CBar.set_label(r'$Log_{10} PSD (s^3/cm^6)$', x=0, y=0.55, rotation=90, fontsize=10)
        #
        Time=Duration[0]
        TimeStr = str(Time[0]) + '-' + str(Time[1]) + '-' + str(Time[2]) + '/' + \
                  str(Time[3]) + ':' + str(Time[4]) + ':' + str(Time[5]) + '.' + \
                  str(Time[6])+ '+' + str(tt-EpochArrToCal[0])
        AX.text(4, 16.5, TimeStr, fontsize=10)
        if not os.path.exists(FigurePath):
            os.makedirs(FigurePath)
        FigureName = FigurePath + 'MMS1_trans' + str(Time[0]) + str(Time[1]) + str(Time[2]) + \
                     '-' + str(Time[3]) + str(Time[4]) + str(Time[5]) + '.' + str(Time[6]) + \
                     '_' +  '+' + str(tt - EpochArrToCal[0])+\
                     '.png'

        plt.savefig(FigureName)
        plt.show()
        plt.close()
    return()
    #return(DistTT20)


def MomentsAll(Duration, FileDist):
    ###get data==========================
    InfoDist=FileDist.cdf_info()
    EpochData=FileDist.varget('Epoch')
    PhiData=FileDist.varget('mms1_des_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_des_phi_delta_brst')
    DistData=FileDist.varget('mms1_des_dist_brst')
    DistErrData=FileDist.varget('mms1_des_disterr_brst')
    ThetaData=FileDist.varget('mms1_des_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_des_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_des_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_des_energy_delta_brst')
    ###===================================get data
    ###convert data to new data in SI units==========================
    DistData=DistData*10**12###dist in s3/m6
    DistErrData=DistErrData*10**12###dist in s3/m6
    EnergyChannel=EnergyChannel*eV2J###energy in J
    EnergyChannelDelta=EnergyChannelDelta*eV2J###energy in J
    ThetaData=ThetaData*math.pi/180###angle in radian
    ThetaDataDelta=ThetaDataDelta*2*math.pi/180###angle in radian
    PhiData=PhiData*math.pi/180###angle in radian
    PhiDataDelta=PhiDataDelta*2*math.pi/180###angle in radian
    ###
    ###Get the epoch slice==================================================
    EpochStartType2000=cdflib.cdfepoch.compute_tt2000(Duration[0])
    EpochEndType2000=cdflib.cdfepoch.compute_tt2000(Duration[1])
    EpochTplStart=np.where(EpochData >= EpochStartType2000)
    EpochTplEnd=np.where(EpochData <= EpochEndType2000)
    EpochArr1=EpochTplStart[0]
    EpochArr2=EpochTplEnd[0]
    EpochArrToCal=np.intersect1d(EpochArr1, EpochArr2)
    ####=================================Get the epoch slice

    den=[]
    vx=[]
    vy=[]
    vz=[]
    vtotal=[]
    pth=[]
    for tt in EpochArrToCal:
        ###
        DistTT = DistData[tt]
        DistErrTT = DistErrData[tt]
        EnergyChannelTT = EnergyChannel[tt]
        EnergyChannelDeltaTT = EnergyChannelDelta[tt]
        ThetaTT = ThetaData
        ThetaDelta = ThetaDataDelta[0]
        PhiTT = PhiData[tt]
        PhiDelta = PhiDataDelta[0]
        ###
        ###
        denTT = 0
        nByVxTT = 0
        nByVyTT = 0
        nByVzTT = 0
        pxxItem1 = 0
        pyyItem1 = 0
        pzzItem1 = 0
        ###
        ###
        # Energy20TT=EnergyChannelTT[:,:,20]
        for ii in np.arange(32):
            EnergyII = EnergyChannelTT[ii]
            EnergyDeltaII = EnergyChannelDeltaTT[ii]
            velII = math.sqrt(2 * EnergyII / mp)
            velIISqrt = velII ** 2
            velIICube = velII ** 3
            velIIQuar = velII ** 4
            EnergyIILeft = EnergyII - EnergyDeltaII
            EnergyIIRight = EnergyII + EnergyDeltaII
            velIILeft = math.sqrt(2 * EnergyIILeft / mp)
            velIIRight = math.sqrt(2 * EnergyIIRight / mp)
            velDeltaII = velIIRight - velIILeft
            for jj in np.arange(32):
                PhiJJ = PhiTT[jj]
                for kk in np.arange(16):
                    ThetaKK = ThetaTT[kk]
                    denTT += ThetaDelta * PhiDelta * velIISqrt * velDeltaII * math.sin(ThetaKK) * DistTT[jj, kk, ii]
                    nByVxTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * (
                        math.sin(ThetaKK)) ** 2 * math.cos(PhiJJ) * DistTT[jj, kk, ii]
                    nByVyTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * (
                        math.sin(ThetaKK)) ** 2 * math.sin(PhiJJ) * DistTT[jj, kk, ii]
                    nByVzTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * math.sin(ThetaKK) * math.cos(
                        ThetaKK) * DistTT[jj, kk, ii]
                    pxxItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * (math.sin(ThetaKK)) ** 3 * (
                        math.cos(PhiJJ)) ** 2 * DistTT[jj, kk, ii]
                    pyyItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * (math.sin(ThetaKK)) ** 3 * (
                        math.sin(PhiJJ)) ** 2 * DistTT[jj, kk, ii]
                    pzzItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * math.sin(ThetaKK) * (
                        math.cos(ThetaKK)) ** 2 * DistTT[jj, kk, ii]

        den = np.append(den, denTT)

        vxTT = nByVxTT / denTT
        vyTT = nByVyTT / denTT
        vzTT = nByVzTT / denTT
        vtotalTT = math.sqrt(vxTT ** 2 + vyTT ** 2 + vzTT ** 2)
        vx=np.append(vx, vxTT)
        vy=np.append(vy, vyTT)
        vz=np.append(vz, vzTT)
        vtotal=np.append(vtotal, vtotalTT)

        pxxTT = pxxItem1 - mp * denTT * vxTT ** 2
        pyyTT = pyyItem1 - mp * denTT * vyTT ** 2
        pzzTT = pzzItem1 - mp * denTT * vzTT ** 2
        pTT = (pxxTT + pyyTT + pzzTT) / 3
        pth=np.append(pth, pTT)
        #
    densityM = den  ###number density in m**(-3)
    densityCM = den * 10 ** (-6)  ###number density in cm*(-3)
    vxMS = vx  ###velocity in m/s
    vyMS = vy  ###velocity in m/s
    vzMS = vz  ###velocity in m/s
    vtotalMS = vtotal  ###velocity in m/s
    vxKMS = vx * 10 ** (-3)  ###velocity in km/s
    vyKMS = vy * 10 ** (-3)  ###velocity in km/s
    vzKMS = vz * 10 ** (-3)  ###velocity in km/s
    vtotalKMS = vtotal * 10 ** (-3)  ###velocity in km/s
    pthPa = pth  ###p in Pa
    pthnPa = pthPa * 10 ** 9  ###p in nPa

    tempK = pthPa / (densityM * kb)
    tempeV = (pthPa / densityM) / eV2J

    EpochDist=EpochData[EpochArrToCal]
    return (densityCM, vxKMS, vyKMS, vzKMS, pthnPa, tempeV, EpochDist)


def MomentsNo(Duration, FileDist, Theta, Phi):
    ###init============================
    PhiNo=np.array(Phi)
    ThetaNo=np.array(Theta)
    ###===================================init
    ###get data==========================
    InfoDist=FileDist.cdf_info()
    EpochData=FileDist.varget('Epoch')
    PhiData=FileDist.varget('mms1_des_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_des_phi_delta_brst')
    DistData=FileDist.varget('mms1_des_dist_brst')
    DistErrData=FileDist.varget('mms1_des_disterr_brst')
    ThetaData=FileDist.varget('mms1_des_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_des_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_des_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_des_energy_delta_brst')
    ###===================================get data
    ###convert data to new data in SI units==========================
    DistData=DistData*10**12###dist in s3/m6
    DistErrData=DistErrData*10**12###dist in s3/m6
    EnergyChannel=EnergyChannel*eV2J###energy in J
    EnergyChannelDelta=EnergyChannelDelta*eV2J###energy in J
    ThetaData=ThetaData*math.pi/180###angle in radian
    ThetaDataDelta=ThetaDataDelta*2*math.pi/180###angle in radian
    PhiData=PhiData*math.pi/180###angle in radian
    PhiDataDelta=PhiDataDelta*2*math.pi/180###angle in radian
    ###
    ###Get the epoch slice==================================================
    EpochStartType2000=cdflib.cdfepoch.compute_tt2000(Duration[0])
    EpochEndType2000=cdflib.cdfepoch.compute_tt2000(Duration[1])
    EpochTplStart=np.where(EpochData >= EpochStartType2000)
    EpochTplEnd=np.where(EpochData <= EpochEndType2000)
    EpochArr1=EpochTplStart[0]
    EpochArr2=EpochTplEnd[0]
    EpochArrToCal=np.intersect1d(EpochArr1, EpochArr2)
    ####=================================Get the epoch slice

    den=[]
    vx=[]
    vy=[]
    vz=[]
    vtotal=[]
    pth=[]
    for tt in EpochArrToCal:
        ###
        DistTT = DistData[tt]
        DistErrTT = DistErrData[tt]
        EnergyChannelTT = EnergyChannel[tt]
        EnergyChannelDeltaTT = EnergyChannelDelta[tt]
        ThetaTT = ThetaData
        ThetaDelta = ThetaDataDelta[0]
        PhiTT = PhiData[tt]
        PhiDelta = PhiDataDelta[0]
        ###
        ###
        denTT = 0
        nByVxTT = 0
        nByVyTT = 0
        nByVzTT = 0
        pxxItem1 = 0
        pyyItem1 = 0
        pzzItem1 = 0
        ###
        ###
        # Energy20TT=EnergyChannelTT[:,:,20]
        for ii in np.arange(32):
            EnergyII = EnergyChannelTT[ii]
            EnergyDeltaII = EnergyChannelDeltaTT[ii]
            velII = math.sqrt(2 * EnergyII / mp)
            velIISqrt = velII ** 2
            velIICube = velII ** 3
            velIIQuar = velII ** 4
            EnergyIILeft = EnergyII - EnergyDeltaII
            EnergyIIRight = EnergyII + EnergyDeltaII
            velIILeft = math.sqrt(2 * EnergyIILeft / mp)
            velIIRight = math.sqrt(2 * EnergyIIRight / mp)
            velDeltaII = velIIRight - velIILeft
            for jj in np.arange(32):
                PhiJJ = PhiTT[jj]
                for kk in np.arange(16):
                    ThetaKK = ThetaTT[kk]
                    ###No===========
                    if jj in PhiNo:
                        if kk in ThetaNo:
                            DistTT[jj, kk, :] = 0
                    ###=============No
                    denTT += ThetaDelta * PhiDelta * velIISqrt * velDeltaII * math.sin(ThetaKK) * DistTT[jj, kk, ii]
                    nByVxTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * (
                        math.sin(ThetaKK)) ** 2 * math.cos(PhiJJ) * DistTT[jj, kk, ii]
                    nByVyTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * (
                        math.sin(ThetaKK)) ** 2 * math.sin(PhiJJ) * DistTT[jj, kk, ii]
                    nByVzTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * math.sin(ThetaKK) * math.cos(
                        ThetaKK) * DistTT[jj, kk, ii]
                    pxxItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * (math.sin(ThetaKK)) ** 3 * (
                        math.cos(PhiJJ)) ** 2 * DistTT[jj, kk, ii]
                    pyyItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * (math.sin(ThetaKK)) ** 3 * (
                        math.sin(PhiJJ)) ** 2 * DistTT[jj, kk, ii]
                    pzzItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * math.sin(ThetaKK) * (
                        math.cos(ThetaKK)) ** 2 * DistTT[jj, kk, ii]

        den = np.append(den, denTT)

        vxTT = nByVxTT / denTT
        vyTT = nByVyTT / denTT
        vzTT = nByVzTT / denTT
        vtotalTT = math.sqrt(vxTT ** 2 + vyTT ** 2 + vzTT ** 2)
        vx=np.append(vx, vxTT)
        vy=np.append(vy, vyTT)
        vz=np.append(vz, vzTT)
        vtotal=np.append(vtotal, vtotalTT)

        pxxTT = pxxItem1 - mp * denTT * vxTT ** 2
        pyyTT = pyyItem1 - mp * denTT * vyTT ** 2
        pzzTT = pzzItem1 - mp * denTT * vzTT ** 2
        pTT = (pxxTT + pyyTT + pzzTT) / 3
        pth=np.append(pth, pTT)
        #
    densityM = den  ###number density in m**(-3)
    densityCM = den * 10 ** (-6)  ###number density in cm*(-3)
    vxMS = vx  ###velocity in m/s
    vyMS = vy  ###velocity in m/s
    vzMS = vz  ###velocity in m/s
    vtotalMS = vtotal  ###velocity in m/s
    vxKMS = vx * 10 ** (-3)  ###velocity in km/s
    vyKMS = vy * 10 ** (-3)  ###velocity in km/s
    vzKMS = vz * 10 ** (-3)  ###velocity in km/s
    vtotalKMS = vtotal * 10 ** (-3)  ###velocity in km/s
    pthPa = pth  ###p in Pa
    pthnPa = pthPa * 10 ** 9  ###p in nPa

    tempK = pthPa / (densityM * kb)
    tempeV = (pthPa / densityM) / eV2J

    return (densityCM, vxKMS, vyKMS, vzKMS, pthnPa, tempeV)


def MomentsYes(Duration, FileDist, Theta, Phi):
    ###init============================
    PhiNo=np.array(Phi)
    ThetaNo=np.array(Theta)
    ###===================================init
    ###get data==========================
    InfoDist=FileDist.cdf_info()
    EpochData=FileDist.varget('Epoch')
    PhiData=FileDist.varget('mms1_des_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_des_phi_delta_brst')
    DistData=FileDist.varget('mms1_des_dist_brst')
    DistErrData=FileDist.varget('mms1_des_disterr_brst')
    ThetaData=FileDist.varget('mms1_des_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_des_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_des_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_des_energy_delta_brst')
    ###===================================get data
    ###convert data to new data in SI units==========================
    DistData=DistData*10**12###dist in s3/m6
    DistErrData=DistErrData*10**12###dist in s3/m6
    EnergyChannel=EnergyChannel*eV2J###energy in J
    EnergyChannelDelta=EnergyChannelDelta*eV2J###energy in J
    ThetaData=ThetaData*math.pi/180###angle in radian
    ThetaDataDelta=ThetaDataDelta*2*math.pi/180###angle in radian
    PhiData=PhiData*math.pi/180###angle in radian
    PhiDataDelta=PhiDataDelta*2*math.pi/180###angle in radian

    ###Get the epoch slice==================================================
    EpochStartType2000=cdflib.cdfepoch.compute_tt2000(Duration[0])
    EpochEndType2000=cdflib.cdfepoch.compute_tt2000(Duration[1])
    EpochTplStart=np.where(EpochData >= EpochStartType2000)
    EpochTplEnd=np.where(EpochData <= EpochEndType2000)
    EpochArr1=EpochTplStart[0]
    EpochArr2=EpochTplEnd[0]
    EpochArrToCal=np.intersect1d(EpochArr1, EpochArr2)
    ####=================================Get the epoch slice
    ###
    den=[]
    vx=[]
    vy=[]
    vz=[]
    vtotal=[]
    pth=[]
    for tt in EpochArrToCal:
        ###
        DistTT = DistData[tt]
        DistErrTT = DistErrData[tt]
        EnergyChannelTT = EnergyChannel[tt]
        EnergyChannelDeltaTT = EnergyChannelDelta[tt]
        ThetaTT = ThetaData
        ThetaDelta = ThetaDataDelta[0]
        PhiTT = PhiData[tt]
        PhiDelta = PhiDataDelta[0]
        ###
        ###
        denTT = 0
        nByVxTT = 0
        nByVyTT = 0
        nByVzTT = 0
        pxxItem1 = 0
        pyyItem1 = 0
        pzzItem1 = 0
        ###
        ###
        for ii in np.arange(32):
            EnergyII = EnergyChannelTT[ii]
            EnergyDeltaII = EnergyChannelDeltaTT[ii]
            velII = math.sqrt(2 * EnergyII / mp)
            velIISqrt = velII ** 2
            velIICube = velII ** 3
            velIIQuar = velII ** 4
            EnergyIILeft = EnergyII - EnergyDeltaII
            EnergyIIRight = EnergyII + EnergyDeltaII
            velIILeft = math.sqrt(2 * EnergyIILeft / mp)
            velIIRight = math.sqrt(2 * EnergyIIRight / mp)
            velDeltaII = velIIRight - velIILeft
            for jj in np.arange(32):
                PhiJJ = PhiTT[jj]
                for kk in np.arange(16):
                    ThetaKK = ThetaTT[kk]
                    ###No===========
                    if jj in PhiNo:
                        if kk in ThetaNo:
                            denTT += ThetaDelta * PhiDelta * velIISqrt * velDeltaII * math.sin(ThetaKK) * DistTT[
                                jj, kk, ii]
                            nByVxTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * (
                                math.sin(ThetaKK)) ** 2 * math.cos(PhiJJ) * DistTT[jj, kk, ii]
                            nByVyTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * (
                                math.sin(ThetaKK)) ** 2 * math.sin(PhiJJ) * DistTT[jj, kk, ii]
                            nByVzTT += (-1) * ThetaDelta * PhiDelta * velIICube * velDeltaII * math.sin(
                                ThetaKK) * math.cos(ThetaKK) * DistTT[jj, kk, ii]
                            pxxItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * (
                                math.sin(ThetaKK)) ** 3 * (math.cos(PhiJJ)) ** 2 * DistTT[jj, kk, ii]
                            pyyItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * (
                                math.sin(ThetaKK)) ** 3 * (math.sin(PhiJJ)) ** 2 * DistTT[jj, kk, ii]
                            pzzItem1 += mp * ThetaDelta * PhiDelta * velIIQuar * velDeltaII * math.sin(ThetaKK) * (
                                math.cos(ThetaKK)) ** 2 * DistTT[jj, kk, ii]

        den=np.append(den, denTT)

        vxTT = nByVxTT / denTT
        vyTT = nByVyTT / denTT
        vzTT = nByVzTT / denTT
        vtotalTT = math.sqrt(vxTT ** 2 + vyTT ** 2 + vzTT ** 2)
        vx=np.append(vx, vxTT)
        vy=np.append(vy, vyTT)
        vz=np.append(vz, vzTT)
        vtotal=np.append(vtotal, vtotalTT)

        pxxTT = pxxItem1 - mp * denTT * vxTT ** 2
        pyyTT = pyyItem1 - mp * denTT * vyTT ** 2
        pzzTT = pzzItem1 - mp * denTT * vzTT ** 2
        pTT = (pxxTT + pyyTT + pzzTT) / 3
        pth=np.append(pth, pTT)
        #
    densityM = den  ###number density in m**(-3)
    densityCM = den * 10 ** (-6)  ###number density in cm*(-3)
    vxMS = vx  ###velocity in m/s
    vyMS = vy  ###velocity in m/s
    vzMS = vz  ###velocity in m/s
    vtotalMS = vtotal  ###velocity in m/s
    vxKMS = vx * 10 ** (-3)  ###velocity in km/s
    vyKMS = vy * 10 ** (-3)  ###velocity in km/s
    vzKMS = vz * 10 ** (-3)  ###velocity in km/s
    vtotalKMS = vtotal * 10 ** (-3)  ###velocity in km/s
    pthPa = pth  ###p in Pa
    pthnPa = pthPa * 10 ** 9  ###p in nPa

    tempK = pthPa / (densityM * kb)
    tempeV = (pthPa / densityM) / eV2J

    return (densityCM, vxKMS, vyKMS, vzKMS, pthnPa, tempeV)




###=======================================================
PhiInput =[]
ThetaInput =[]

###test
TimeInputStart = [2017, 12, 1, 14, 41, 5, 0]
TimeInputEnd = [2017, 12, 1, 14, 41, 5, 200]
DurationInput = [TimeInputStart, TimeInputEnd]

FileInput = \
    cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/des-dist/2017/12/01/mms1_fpi_brst_l2_des-dist_20171201143933_v3.3.0.cdf')
FigurePathInput = 'D:/OneDrive - mail.sdu.edu.cn/ResearchWork/code&figure/'
#
SkyMap=SkyMap(DurationInput,FileInput,ThetaInput,PhiInput,FigurePathInput)

EnergyChannel=FileInput.varget('mms1_des_energy_brst')
EE=EnergyChannel[0]
np.set_printoptions(precision=3, suppress=True)
print('\n {}'.format(EE))


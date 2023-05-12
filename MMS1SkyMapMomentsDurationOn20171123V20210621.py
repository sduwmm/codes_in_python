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
'mms1_dis_errorflags_brst', 'mms1_dis_compressionloss_brst', 
'mms1_dis_steptable_parity_brst', 'mms1_dis_startdelphi_count_brst', 
'mms1_dis_startdelphi_angle_brst', 
'mms1_dis_phi_brst', 'mms1_dis_phi_delta_brst', 
'mms1_dis_dist_brst', 'mms1_dis_disterr_brst', 
'mms1_dis_sector_despinp_brst', 'mms1_dis_sector_label_brst', 
'mms1_dis_pixel_label_brst', 'mms1_dis_energy_label_brst', 
'mms1_dis_theta_brst', 'mms1_dis_theta_delta_brst', 
'mms1_dis_energy_brst', 'mms1_dis_energy_delta_brst']
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
    PhiData=FileDist.varget('mms1_dis_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_dis_phi_delta_brst')
    DistData=FileDist.varget('mms1_dis_dist_brst')
    DistErrData=FileDist.varget('mms1_dis_disterr_brst')
    ThetaData=FileDist.varget('mms1_dis_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_dis_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_dis_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_dis_energy_delta_brst')
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
        for ee in np.arange(0, 32, 1):  # without bins below 5eV 2021-06-17

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
        CMap = AX.pcolor(DistSumLg, cmap='jet', vmin=-26, vmax=-20)
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
        #plt.show()
        #plt.close()
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
    PhiData=FileDist.varget('mms1_dis_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_dis_phi_delta_brst')
    DistData=FileDist.varget('mms1_dis_dist_brst')
    DistErrData=FileDist.varget('mms1_dis_disterr_brst')
    ThetaData=FileDist.varget('mms1_dis_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_dis_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_dis_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_dis_energy_delta_brst')
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
        for ee in np.arange(0, 32, 1):  # without bins below 5eV 2021-06-17

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
    PhiData=FileDist.varget('mms1_dis_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_dis_phi_delta_brst')
    DistData=FileDist.varget('mms1_dis_dist_brst')
    DistErrData=FileDist.varget('mms1_dis_disterr_brst')
    ThetaData=FileDist.varget('mms1_dis_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_dis_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_dis_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_dis_energy_delta_brst')
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
    PhiData=FileDist.varget('mms1_dis_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_dis_phi_delta_brst')
    DistData=FileDist.varget('mms1_dis_dist_brst')
    DistErrData=FileDist.varget('mms1_dis_disterr_brst')
    ThetaData=FileDist.varget('mms1_dis_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_dis_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_dis_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_dis_energy_delta_brst')
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
    PhiData=FileDist.varget('mms1_dis_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_dis_phi_delta_brst')
    DistData=FileDist.varget('mms1_dis_dist_brst')
    DistErrData=FileDist.varget('mms1_dis_disterr_brst')
    ThetaData=FileDist.varget('mms1_dis_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_dis_theta_delta_brst')
    EnergyChannel=FileDist.varget('mms1_dis_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_dis_energy_delta_brst')
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
PhiInput =[0,1,30,31]
ThetaInput =np.arange(4, 11, 1)
#momFore = MomentsNo(TimeInput, FileInput, PhiInput, ThetaInput)
###test
TimeInputStart = [2017, 11, 23, 4, 53, 20, 0]
TimeInputEnd = [2017, 11, 23, 4, 54, 20, 0]
DurationInput = [TimeInputStart, TimeInputEnd]

FileInput = \
    cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2017/11/23/mms1_fpi_brst_l2_dis-dist_20171123045143_v3.3.0.cdf')
FigurePathInput = 'D:/OneDrive - mail.sdu.edu.cn/ResearchWork/code&figure/figure2021/skyMap/'
#
#SkyMap=SkyMap(DurationInput,FileInput,ThetaInput,PhiInput,FigurePathInput)
#SkyMapTrans=SkyMapTrans(DurationInput,FileInput,ThetaInput,PhiInput,FigurePathInput)
momAll=MomentsAll(DurationInput,FileInput)
momNo=MomentsNo(DurationInput,FileInput,ThetaInput,PhiInput)
momYes=MomentsYes(DurationInput,FileInput,ThetaInput,PhiInput)
#
DenAll = momAll[0]
EpochDist=momAll[6]
DenNo = momNo[0]####
DenYes = momYes[0]#####SW beam
PthAll = momAll[4]
PthNo = momNo[4]
PthYes = momYes[4]

#####==================================================================================================================
#####==================================================================================================================
FileMom = \
    cdflib.CDF('D:\\data\\mms\\mms1\\fpi\\brst\\l2\\dis-moms\\2017\\11\\23\mms1_fpi_brst_l2_dis-moms_20171123045143_v3.3.0.cdf')
InfoMom = FileMom.cdf_info()
EpochMom = FileMom.varget('Epoch')
DenMom = FileMom.varget('mms1_dis_numberdensity_brst')
DenErrMom = FileMom.varget('mms1_dis_numberdensity_err_brst')
#BulkvDBCSMom = FileMom.varget('mms1_dis_bulkv_dbcs_brst')####'mms1_dis_bulkv_spintone_dbcs_brst'
TempParaMom = FileMom.varget('mms1_dis_temppara_brst')
TempPerpMom = FileMom.varget('mms1_dis_tempperp_brst')
EnerSpecMom = FileMom.varget('mms1_dis_energyspectr_omni_brst')

###Get the epoch slice==================================================
EpochStartType2000 = cdflib.cdfepoch.compute_tt2000(DurationInput[0])
EpochEndType2000 = cdflib.cdfepoch.compute_tt2000(DurationInput[1])
EpochTplStartMom = np.where(EpochMom >= EpochStartType2000)
EpochTplEndMom = np.where(EpochMom <= EpochEndType2000)
EpochArr1Mom = EpochTplStartMom[0]
EpochArr2Mom = EpochTplEndMom[0]
EpochArrToCalMom = np.intersect1d(EpochArr1Mom, EpochArr2Mom)
####=================================Get the epoch slice
#
EpochMomToPlot = EpochMom[EpochArrToCalMom]
DenOff = DenMom[EpochArrToCalMom]
#BulkvDBCSOff = BulkvDBCSMom[EpochArrToCalMom]
TempParaOff = TempParaMom[EpochArrToCalMom]
TempPerpOff = TempPerpMom[EpochArrToCalMom]
#
#vxOff = BulkvDBCSOff[:,0]
#vyOff = BulkvDBCSOff[:,1]
#vzOff = BulkvDBCSOff[:,2]
#vtotalOff = np.power((np.power(vxOff,2) + np.power(vyOff,2) + np.power(vzOff,2)), 1/2)
TempAvgOff = 1/3*TempParaOff + 2/3*TempPerpOff
TempAvgOffJ = TempAvgOff * eV2J
DenOffM = DenOff*10**6#in m**(-6)
PthOffPa = DenOffM * TempAvgOffJ
PthOffnPa = PthOffPa*10**9
###
FileFgm = \
    cdflib.CDF('D:\\data\\mms\\mms1\\fgm\\brst\\l2\\2017\\11\\23\\mms1_fgm_brst_l2_20171123045143_v5.113.0.cdf')
InfoFgm = FileFgm.cdf_info()
EpochFgm = FileFgm.varget('Epoch')
BgseFgm = FileFgm.varget('mms1_fgm_b_gse_brst_l2')
###Get the epoch slice==================================================
EpochStartType2000 = cdflib.cdfepoch.compute_tt2000(DurationInput[0])
EpochEndType2000 = cdflib.cdfepoch.compute_tt2000(DurationInput[1])
EpochTplStartFgm = np.where(EpochFgm >= EpochStartType2000)
EpochTplEndFgm = np.where(EpochFgm <= EpochEndType2000)
EpochArr1Fgm = EpochTplStartFgm[0]
EpochArr2Fgm = EpochTplEndFgm[0]
EpochArrToCalFgm = np.intersect1d(EpochArr1Fgm, EpochArr2Fgm)
####=================================Get the epoch slice
EpochFgmToPlot = EpochFgm[EpochArrToCalFgm]
Bgse = BgseFgm[EpochArrToCalFgm]
###===============================================================================
Fig = plt.figure(figsize=(10,10), dpi=300)
Fig.subplots_adjust(top=0.90)
Fig.subplots_adjust(bottom=0.07)
Fig.subplots_adjust(left=0.13)
Fig.subplots_adjust(right=0.93)
ax1=plt.subplot(5,1,1)
ax2=plt.subplot(5,1,2)
ax3=plt.subplot(5,1,3)
ax4=plt.subplot(5,1,4)
ax5=plt.subplot(5,1,5)


xlimFgm=[min(EpochFgmToPlot), max(EpochFgmToPlot)]
ax1.plot(EpochFgmToPlot, Bgse[:,0], 'blue', label = 'Bx')
ax1.plot(EpochFgmToPlot, Bgse[:,1], 'green', label = 'By')
ax1.plot(EpochFgmToPlot, Bgse[:,2], 'red',label='Bz')
ax1.plot(EpochFgmToPlot, Bgse[:,3], 'black',label='|B|')
ax1.legend(loc = 1, fontsize = 12, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax1.set_xlim(xlimFgm)
ax1.set_title('MMS1 2017-11-23',fontsize=20)
ax1.axes.get_xaxis().set_visible(False)
ax1.set_ylabel('B FGM [nT]', fontsize=16)

xlimDist=[min(EpochDist),max(EpochDist)]
ax2.plot(EpochDist,DenAll, color='green',linestyle = '-',label='Nions-1beam')
ax2.plot(EpochMomToPlot, DenOff, color = 'black',label='Nions-official')
ax2.legend(loc = 1, fontsize = 12, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax2.set_xlim(xlimDist)
ax2.set_xticks([])
ax2.set_ylabel('Nions', fontsize=16)

ax3.plot(EpochDist,DenAll, color = 'green',label='Nions-1beam')
DenSum = DenNo+DenYes
ax3.plot(EpochDist, DenSum,color='blue',label='Nions-2beams')
ax3.legend(loc = 1, fontsize = 12, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax3.set_xlim(xlimDist)
ax3.set_xticks([])
ax3.set_ylabel('Nions', fontsize=16)

ax4.plot(EpochDist,PthAll, color = 'orange', label='Pth-1beam')
ax4.plot(EpochMomToPlot,PthOffnPa, color = 'pink',label = 'Pth-official')
ax4.legend(loc = 1, fontsize = 12, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax4.set_xlim(xlimDist)
ax4.set_ylim([0.05,0.3])
ax4.set_xticks([])
ax4.set_ylabel('pth [nPa]', fontsize=16)

PthSum=PthNo+PthYes
ax5.plot(EpochDist, PthSum, color = 'deeppink', label = 'Pth-2beams')
ax5.plot(EpochMomToPlot, PthOffnPa, color = 'pink',label='Pth-official')
ax5.legend(loc = 1, fontsize = 12, markerscale = 1.0, handlelength = 2, handleheight = 0.5)
ax5.set_xlim(xlimDist)
ax5.set_ylim([0.05,0.3])
ax5.set_xticks([])
ax5.set_ylabel('pth [nPa]', fontsize=16)
ax5.tick_params(labelsize=15)
#ax5.text(EpochDist[len(EpochDist)-10], np.median(PthSum), '(e)', fontsize = 22)

xtick1 = [2017, 11, 23, 4, 53, 30, 000]
xtick2 = [2017, 11, 23, 4, 54, 0, 000]
xtick1 = cdflib.cdfepoch.compute_tt2000(xtick1)
xtick2 = cdflib.cdfepoch.compute_tt2000(xtick2)
epochTick = [xtick1, xtick2]
strTick = ['04:53:30', '04:54:00']
ax5.set_xticks(epochTick)
ax5.set_xticklabels(strTick, fontsize = 15)


figurePath = 'D:/OneDrive - mail.sdu.edu.cn/'
if not os.path.exists(figurePath):
    os.makedirs(figurePath)

plt.savefig(figurePath+'/beam-identify-recalculate.png',dpi = 200)
plt.show()
#plt.close()

pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))
###get data from dis-mom


EnergyChannel=FileInput.varget('mms1_dis_energy_brst')
EE=EnergyChannel[0]
np.set_printoptions(precision=3, suppress=True)
print('\n {}'.format(EE))

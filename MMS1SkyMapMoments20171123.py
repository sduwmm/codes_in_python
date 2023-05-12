'''
Created on Tue Apr 9 11:09:20 2019
Changed on Tue Mar 23 17:33:39 2021
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

import cdflib
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import time

pyStartTime = time.time()
def SkyMap(Time, FileDist, Theta, Phi, FigurePath):
    ###init============================
    ThetaNo=np.array(Theta)
    PhiNo = np.array(Phi)
    ###===================================init
    ###get data==========================
    InfoDist=FileDist.cdf_info()
    EpochData=FileDist.varget('Epoch')
    DistData=FileDist.varget('mms1_dis_dist_brst')
    DistErrData=FileDist.varget('mms1_dis_disterr_brst')
    EnergyChannel=FileDist.varget('mms1_dis_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_dis_energy_delta_brst')
    ThetaData=FileDist.varget('mms1_dis_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_dis_theta_delta_brst')
    PhiData=FileDist.varget('mms1_dis_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_dis_phi_delta_brst')
    ###===================================get data
    ###
    EpochToCalculateType2000=cdflib.cdfepoch.compute_tt2000(Time)
    EpochIndTpl=np.where(EpochData>=EpochToCalculateType2000)
    EpochIndArr=EpochIndTpl[0]
    EpochInd=EpochIndArr[0]
    ###
    DistTT=DistData[EpochInd]
    DistErrTT=DistErrData[EpochInd]
    EnergyChannelTT=EnergyChannel[EpochInd]
    EnergyChannelDeltaTT=EnergyChannelDelta[EpochInd]
    #ThetaTT=ThetaData
    #PhiTT=PhiData[EpochInd]
    #
    DistAlongEnergySum=np.zeros([16,32])
    ###
    for ee in np.arange(0,32,1):
        EnergyInd=ee
        DistTTEE=DistTT[ee]#take out the dist of certain energy channel
        DistErrTTEE=DistErrTT[ee]#take out the dist error of certain energy channel
        #
        #'''
        for jj in ThetaNo:
            for kk in PhiNo:
                if jj in ThetaNo:
                    if kk in PhiNo:
                        DistTTEE[jj, kk] = 0
        #'''
        #
        #
        #DistTTEE[DistTTEE<=DistErrTTEE]=0
        #
        DistAlongEnergySum=DistAlongEnergySum+DistTTEE
        #
        DistSumLg=np.log10(DistAlongEnergySum)
    #
    Figure=plt.figure(figsize=(5,4),dpi=100)
    AX=plt.subplot(1,1,1)
    CMap=AX.pcolor(DistSumLg,cmap='jet',vmin=-26,vmax=-20)
    CBar=Figure.colorbar(CMap,ax=AX)
    CBar.set_label(r'$Log_{10} PSD (s^3/cm^6)$', x=0, y=0.55, rotation=90, fontsize=10)
    #
    TimeStr = str(Time[0]) + '-' + str(Time[1]) + '-' + str(Time[2]) + '/' + \
              str(Time[3]) + ':' + str(Time[4]) + ':' + str(Time[5]) + '.' + \
              str(Time[6])
    AX.text(4, 16.5, TimeStr, fontsize=10)
    if not os.path.exists(FigurePath):
        os.makedirs(FigurePath)
    FigureName = FigurePath + 'MMS1_' + str(Time[0]) + str(Time[1]) + str(Time[2]) + \
                 '-' + str(Time[3]) + str(Time[4]) + str(Time[5]) + '.' + str(Time[6]) + \
                 '.png'

    plt.savefig(FigureName)
    plt.show()
    #plt.close()
    return(InfoDist, DistData)


    #

def MomentsNo(Time, FileDist, Theta, Phi):
    ###constants===============================
    mp=1.67262158*10**(-27)
    eV2J=1.602176462*10**(-19)
    kb=1.3806503*10**(-23)
    mu0=4*math.pi*10**(-7)
    ###=================================constants
    ###init============================
    PhiNo=np.array(Phi)
    ThetaNo=np.array(Theta)
    ###===================================init
    ###get data==========================
    InfoDist=FileDist.cdf_info()
    EpochData=FileDist.varget('Epoch')
    DistData=FileDist.varget('mms1_dis_dist_brst')
    DistErrData=FileDist.varget('mms1_dis_disterr_brst')
    EnergyChannel=FileDist.varget('mms1_dis_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_dis_energy_delta_brst')
    ThetaData=FileDist.varget('mms1_dis_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_dis_theta_delta_brst')
    PhiData=FileDist.varget('mms1_dis_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_dis_phi_delta_brst')
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
    EpochToCalculateType2000=cdflib.cdfepoch.compute_tt2000(Time)
    EpochIndTpl=np.where(EpochData>=EpochToCalculateType2000)
    EpochIndArr=EpochIndTpl[0]
    EpochInd=EpochIndArr[0]
    ###
    DistTT=DistData[EpochInd]
    DistErrTT=DistErrData[EpochInd]
    EnergyChannelTT=EnergyChannel[EpochInd]
    EnergyChannelDeltaTT=EnergyChannelDelta[EpochInd]
    ThetaTT=ThetaData
    ThetaDelta=ThetaDataDelta[0]
    PhiTT=PhiData[EpochInd]
    PhiDelta=PhiDataDelta[0]
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
        EnergyII=EnergyChannelTT[ii]
        EnergyDeltaII=EnergyChannelDeltaTT[ii]
        velII=math.sqrt(2*EnergyII/mp)
        velIISqrt=velII**2
        velIICube=velII**3
        velIIQuar=velII**4
        EnergyIILeft=EnergyII-EnergyDeltaII
        EnergyIIRight=EnergyII+EnergyDeltaII
        velIILeft=math.sqrt(2*EnergyIILeft/mp)
        velIIRight=math.sqrt(2*EnergyIIRight/mp)
        velDeltaII=velIIRight-velIILeft
        for jj in np.arange(16):
            ThetaJJ=ThetaTT[jj]
            for kk in np.arange(32):
                PhiKK=PhiTT[kk]
                ###No===========
                if jj in ThetaNo:
                    if kk in PhiNo:
                        DistTT[:,jj,kk]=0
                ###=============No
                denTT+=ThetaDelta*PhiDelta*velIISqrt*velDeltaII*math.sin(ThetaJJ)*DistTT[ii,jj,kk]
                nByVxTT+=(-1)*ThetaDelta*PhiDelta*velIICube*velDeltaII*(math.sin(ThetaJJ))**2*math.cos(PhiKK)*DistTT[ii,jj,kk]
                nByVyTT+=(-1)*ThetaDelta*PhiDelta*velIICube*velDeltaII*(math.sin(ThetaJJ))**2*math.sin(PhiKK)*DistTT[ii,jj,kk]
                nByVzTT+=(-1)*ThetaDelta*PhiDelta*velIICube*velDeltaII*math.sin(ThetaJJ)*math.cos(ThetaJJ)*DistTT[ii,jj,kk]
                pxxItem1+=mp*ThetaDelta*PhiDelta*velIIQuar*velDeltaII*(math.sin(ThetaJJ))**3*(math.cos(PhiKK))**2*DistTT[ii,jj,kk]
                pyyItem1+=mp*ThetaDelta*PhiDelta*velIIQuar*velDeltaII*(math.sin(ThetaJJ))**3*(math.sin(PhiKK))**2*DistTT[ii,jj,kk]
                pzzItem1+=mp*ThetaDelta*PhiDelta*velIIQuar*velDeltaII*math.sin(ThetaJJ)*(math.cos(ThetaJJ))**2*DistTT[ii,jj,kk]
    vxTT=nByVxTT/denTT
    vyTT=nByVyTT/denTT
    vzTT=nByVzTT/denTT
    vtotalTT=math.sqrt(vxTT**2+vyTT**2+vzTT**2)
    pxxTT=pxxItem1-mp*denTT*vxTT**2
    pyyTT=pyyItem1-mp*denTT*vyTT**2
    pzzTT=pzzItem1-mp*denTT*vzTT**2
    pTT=(pxxTT+pyyTT+pzzTT)/3
    #
    densityM=denTT  ###number density in m**(-3)
    densityCM=denTT * 10 ** (-6)  ###number density in cm*(-3)
    vxMS=vxTT  ###velocity in m/s
    vyMS=vyTT  ###velocity in m/s
    vzMS=vzTT  ###velocity in m/s
    vtotalMS=vtotalTT  ###velocity in m/s
    vxKMS=vxTT*10**(-3)  ###velocity in km/s
    vyKMS=vyTT*10**(-3)  ###velocity in km/s
    vzKMS=vzTT *10**(-3)  ###velocity in km/s
    vtotalKMS=vtotalTT*10**(-3)  ###velocity in km/s
    pthPa=pTT  ###p in Pa
    pthnPa=pthPa*10**9  ###p in nPa

    tempK=pthPa/(densityM*kb)
    tempeV=(pthPa/densityM)/eV2J

    return (densityCM, vxKMS, vyKMS, vzKMS, pthnPa, tempeV)

def MomentsYes(Time, FileDist, Theta, Phi):
    ###constants===============================
    mp=1.67262158*10**(-27)
    eV2J=1.602176462*10**(-19)
    kb=1.3806503*10**(-23)
    mu0=4*math.pi*10**(-7)
    ###=================================constants
    ###init============================
    PhiNo=np.array(Phi)
    ThetaNo=np.array(Theta)
    ###===================================init
    ###get data==========================
    InfoDist=FileDist.cdf_info()
    EpochData=FileDist.varget('Epoch')
    DistData=FileDist.varget('mms1_dis_dist_brst')
    DistErrData=FileDist.varget('mms1_dis_disterr_brst')
    EnergyChannel=FileDist.varget('mms1_dis_energy_brst')
    EnergyChannelDelta=FileDist.varget('mms1_dis_energy_delta_brst')
    ThetaData=FileDist.varget('mms1_dis_theta_brst')
    ThetaDataDelta=FileDist.varget('mms1_dis_theta_delta_brst')
    PhiData=FileDist.varget('mms1_dis_phi_brst')
    PhiDataDelta=FileDist.varget('mms1_dis_phi_delta_brst')
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
    EpochToCalculateType2000=cdflib.cdfepoch.compute_tt2000(Time)
    EpochIndTpl=np.where(EpochData>=EpochToCalculateType2000)
    EpochIndArr=EpochIndTpl[0]
    EpochInd=EpochIndArr[0]
    ###
    DistTT=DistData[EpochInd]
    DistErrTT=DistErrData[EpochInd]
    EnergyChannelTT=EnergyChannel[EpochInd]
    EnergyChannelDeltaTT=EnergyChannelDelta[EpochInd]
    ThetaTT=ThetaData
    ThetaDelta=ThetaDataDelta[0]
    PhiTT=PhiData[EpochInd]
    PhiDelta=PhiDataDelta[0]
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
        EnergyII=EnergyChannelTT[ii]
        EnergyDeltaII=EnergyChannelDeltaTT[ii]
        velII=math.sqrt(2*EnergyII/mp)
        velIISqrt=velII**2
        velIICube=velII**3
        velIIQuar=velII**4
        EnergyIILeft=EnergyII-EnergyDeltaII
        EnergyIIRight=EnergyII+EnergyDeltaII
        velIILeft=math.sqrt(2*EnergyIILeft/mp)
        velIIRight=math.sqrt(2*EnergyIIRight/mp)
        velDeltaII=velIIRight-velIILeft
        for jj in np.arange(16):
            ThetaJJ=ThetaTT[jj]
            for kk in np.arange(32):
                PhiKK=PhiTT[kk]
                ###No===========
                if jj in ThetaNo:
                    if kk in PhiNo:
                        denTT+=ThetaDelta*PhiDelta*velIISqrt*velDeltaII*math.sin(ThetaJJ)*DistTT[ii,jj,kk]
                        nByVxTT+=(-1)*ThetaDelta*PhiDelta*velIICube*velDeltaII*(math.sin(ThetaJJ))**2*math.cos(PhiKK)*DistTT[ii,jj,kk]
                        nByVyTT+=(-1)*ThetaDelta*PhiDelta*velIICube*velDeltaII*(math.sin(ThetaJJ))**2*math.sin(PhiKK)*DistTT[ii,jj,kk]
                        nByVzTT+=(-1)*ThetaDelta*PhiDelta*velIICube*velDeltaII*math.sin(ThetaJJ)*math.cos(ThetaJJ)*DistTT[ii,jj,kk]
                        pxxItem1+=mp*ThetaDelta*PhiDelta*velIIQuar*velDeltaII*(math.sin(ThetaJJ))**3*(math.cos(PhiKK))**2*DistTT[ii,jj,kk]
                        pyyItem1+=mp*ThetaDelta*PhiDelta*velIIQuar*velDeltaII*(math.sin(ThetaJJ))**3*(math.sin(PhiKK))**2*DistTT[ii,jj,kk]
                        pzzItem1+=mp*ThetaDelta*PhiDelta*velIIQuar*velDeltaII*math.sin(ThetaJJ)*(math.cos(ThetaJJ))**2*DistTT[ii,jj,kk]
    vxTT=nByVxTT/denTT
    vyTT=nByVyTT/denTT
    vzTT=nByVzTT/denTT
    vtotalTT=math.sqrt(vxTT**2+vyTT**2+vzTT**2)
    pxxTT=pxxItem1-mp*denTT*vxTT**2
    pyyTT=pyyItem1-mp*denTT*vyTT**2
    pzzTT=pzzItem1-mp*denTT*vzTT**2
    pTT=(pxxTT+pyyTT+pzzTT)/3
    #
    densityM=denTT  ###number density in m**(-3)
    densityCM=denTT * 10 ** (-6)  ###number density in cm*(-3)
    vxMS=vxTT  ###velocity in m/s
    vyMS=vyTT  ###velocity in m/s
    vzMS=vzTT  ###velocity in m/s
    vtotalMS=vtotalTT  ###velocity in m/s
    vxKMS=vxTT*10**(-3)  ###velocity in km/s
    vyKMS=vyTT*10**(-3)  ###velocity in km/s
    vzKMS=vzTT *10**(-3)  ###velocity in km/s
    vtotalKMS=vtotalTT*10**(-3)  ###velocity in km/s
    pthPa=pTT  ###p in Pa
    pthnPa=pthPa*10**9  ###p in nPa

    tempK=pthPa/(densityM*kb)
    tempeV=(pthPa/densityM)/eV2J

    return (densityCM, vxKMS, vyKMS, vzKMS, pthnPa, tempeV)




###=======================================================
PhiInput =[]# [0,1,2,30,31]#np.arange(0, 3, 1)
ThetaInput =[]# np.arange(3, 12, 1)
#momFore = MomentsNo(TimeInput, FileInput, PhiInput, ThetaInput)
###test
TimeInput = [2017, 11, 23, 4, 53, 55, 0]

FileInput = \
    cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2017/11/23/mms1_fpi_brst_l2_dis-dist_20171123045143_v3.3.0.cdf')
FigurePathInput = 'D:/OneDrive - mail.sdu.edu.cn/work/code&figure/figure2021/skyMap/'

SkyMap=SkyMap(TimeInput,FileInput,ThetaInput,PhiInput,FigurePathInput)
momNo=MomentsNo(TimeInput,FileInput,ThetaInput,PhiInput)
momYes=MomentsYes(TimeInput,FileInput,ThetaInput,PhiInput)


pyEndTime = time.time()
print('==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))
###get data from dis-mom

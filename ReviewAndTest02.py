# -*- coding: utf-8 -*-
"""
Created on Tue Apr 9 11:09:20 2019

@author: Mengmeng Wang at sdu
"""

'''
calculate the zero-, first, second- moments 
of MMS FPI ion distribution data
from the phase space density
in given time,
phi angle and theta angle range,
energy range

reflected beam
for edge

Version 0 
'''
'''
!!!Note
dist: time*energy*theta*phi=epoch_number*32*16*32
'''

import cdflib
import numpy as np
import math
import pyquaternion
import matplotlib.pyplot as plt
##
import time

pyStartTime = time.time()
###
###constants==================
mp   = 1.67262158*10**(-27)
eV2J = 1.602176462*10**(-19)
kb   = 1.3806503*10**(-23)
mu0  = 4*math.pi*10**(-7)
###==================constants
timeInput  = [2017, 11, 16, 12, 13, 40, 500]
normal     = [-0.45, 0.37, -0.80]

loopEnerSW  = np.arange(16,26)
loopPhiSW   = [0, 1, 2, 28, 29, 30, 31]
loopPhiSW   = np.array(loopPhiSW)
loopThetaSW = np.arange(5,14)


fileNameDist   = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2017/11/16/mms1_fpi_brst_l2_dis-dist_20171116121223_v3.3.0.cdf')
infoDist       = fileNameDist.cdf_info()
epochDist      = fileNameDist.varget('epoch')
phiRaw         = fileNameDist.varget('mms1_dis_phi_brst')
phiDeltaRaw    = fileNameDist.varget('mms1_dis_phi_delta_brst')
thetaRaw       = fileNameDist.varget('mms1_dis_theta_brst')
thetaDeltaRaw  = fileNameDist.varget('mms1_dis_theta_delta_brst')
distRaw        = fileNameDist.varget('mms1_dis_dist_brst')
energyRaw      = fileNameDist.varget('mms1_dis_energy_brst')
energyDeltaRaw = fileNameDist.varget('mms1_dis_energy_delta_brst')

epochToCalculateType2000 = cdflib.cdfepoch.compute_tt2000(timeInput)
timeDistIndexTuple           = np.where(epochDist > epochToCalculateType2000)
timeDistIndexArr             = timeDistIndexTuple[0]
timeDistIndexToCal           = timeDistIndexArr[0]


epochDistToCal = epochDist[timeDistIndexToCal]
###convert data to new data in SI units==================================
phi        = phiRaw*math.pi/180###angle in radian
theta      = thetaRaw*math.pi/180###angle in radian
deltaPhi   = phiDeltaRaw*math.pi/180*2###angle in radian
deltaTheta = thetaDeltaRaw*math.pi/180*2###angle in radian
dist       = distRaw*10**12###dist in s3/m6
energy     = energyRaw*eV2J###energy in J
deltaEner  = energyDeltaRaw*eV2J###delta energy in J
###==================================convert data to new data in SI units


phiTT       = phi[timeDistIndexToCal]
thetaTT     = theta
deltaPhi    = deltaPhi[0]
deltaTheta  = deltaTheta[0]
distTT      = dist[timeDistIndexToCal]
energyTT    = energy[timeDistIndexToCal]
deltaEnerTT = deltaEner[timeDistIndexToCal]

loopEner  = np.arange(0, energyTT.size)
loopPhi   = np.arange(0, phiTT.size)
loopTheta = np.arange(0, thetaTT.size)


loopEnerSW  = loopEnerSW
loopPhiSW   = loopPhiSW
loopThetaSW = loopThetaSW

denTT    = 0
nByVxTT  = 0
nByVyTT  = 0
nByVzTT  = 0
pxxItem1 = 0
pyyItem1 = 0
pzzItem1 = 0

for ii in loopEnerSW:
    for jj in loopThetaSW:
        for kk in loopPhiSW:
            distTT[ii, jj, kk] = 0

#for jj in loopTheta:
#    for kk in loopPhi:
#        distTT[4,jj, kk] = 0

for ii in loopEner:
    energyII    = energyTT[ii]
    deltaEnerII = deltaEnerTT[ii]
    velII       = math.sqrt(2*energyII/mp)
    velIISqrt   = velII**2
    velIICube   = velII**3
    velIIQuar   = velII**4
    enerIILeft  = energyII - deltaEnerII
    enerIIRight = energyII + deltaEnerII
    velIILeft   = math.sqrt(2*enerIILeft/mp)
    velIIRight  = math.sqrt(2*enerIIRight/mp)
    velDeltaII  = velIIRight - velIILeft
    for jj in loopTheta:
        thetaJJ = thetaTT[jj]
        for kk in loopPhi:
            phiKK     = phiTT[kk]
            denTT    +=      deltaTheta*deltaPhi*velIISqrt*velDeltaII*math.sin(thetaJJ)*distTT[ii,jj,kk]
            nByVxTT  += (-1)*deltaTheta*deltaPhi*velIICube*velDeltaII*(math.sin(thetaJJ))**2*math.cos(phiKK)*distTT[ii, jj, kk]
            nByVyTT  += (-1)*deltaTheta*deltaPhi*velIICube*velDeltaII*(math.sin(thetaJJ))**2*math.sin(phiKK)*distTT[ii, jj, kk]
            nByVzTT  += (-1)*deltaTheta*deltaPhi*velIICube*velDeltaII*math.sin(thetaJJ)*math.cos(thetaJJ)*distTT[ii, jj, kk]
            pxxItem1 += mp  *deltaTheta*deltaPhi*velIIQuar*velDeltaII*(math.sin(thetaJJ))**3*(math.cos(phiKK))**2*distTT[ii, jj, kk]
            pyyItem1 += mp  *deltaTheta*deltaPhi*velIIQuar*velDeltaII*(math.sin(thetaJJ))**3*(math.sin(phiKK))**2*distTT[ii, jj, kk]
            pzzItem1 += mp  *deltaTheta*deltaPhi*velIIQuar*velDeltaII*math.sin(thetaJJ)*(math.cos(thetaJJ))**2*distTT[ii, jj, kk]

vxTT     = nByVxTT/denTT
vyTT     = nByVyTT/denTT
vzTT     = nByVzTT/denTT
vtotalTT = math.sqrt(vxTT**2 + vyTT**2 + vzTT**2)

pxxTT = pxxItem1 - mp*denTT*vxTT**2
pyyTT = pyyItem1 - mp*denTT*vyTT**2
pzzTT = pzzItem1 - mp*denTT*vzTT**2
pTT   = (pxxTT + pyyTT + pzzTT)/3

densityM  = denTT###number density in m**(-3)
densityCM = denTT*10**(-6)###number density in cm*(-3)
vxMS      = vxTT###velocity in m/s
vyMS      = vyTT###velocity in m/s
vzMS      = vzTT###velocity in m/s
vtotalMS  = vtotalTT###velocity in m/s
vxKMS     = vxTT*10**(-3)###velocity in km/s
vyKMS     = vyTT*10**(-3)###velocity in km/s
vzKMS     = vzTT*10**(-3)###velocity in km/s
vtotalKMS = vtotalTT*10**(-3)###velocity in km/s
pthPa     = pTT###p in Pa
pthnPa    = pthPa*10**9###p in nPa



'''
part 2 convert velocity in DBCS into GSE================================================================================
'''
fileNameMec     = cdflib.CDF('D:/data/mms/mms1/mec/brst/l2/ephts04d/2017/11/16/mms1_mec_brst_l2_ephts04d_20171116121225_v2.2.0.cdf')
infoMec         = fileNameMec.cdf_info()
epochMec        = fileNameMec.varget('Epoch')
quatECI2DBCSArr = fileNameMec.varget('mms1_mec_quat_eci_to_dbcs')
quatECI2GSEArr  = fileNameMec.varget('mms1_mec_quat_eci_to_gse')

epochMecToCalType2000 = cdflib.cdfepoch.compute_tt2000(timeInput)
timeIndexMecTuple  = np.where(epochMec > epochMecToCalType2000)
timeIndexMecArr    = timeIndexMecTuple[0]
timeIndexMecToCal  = timeIndexMecArr[0]

epochMecToCal = epochMec[timeIndexMecToCal]

epochDiff = abs(epochMecToCal - epochDistToCal)
if epochDiff < 150*10**6:
    quatECI2DBCSToCal = quatECI2DBCSArr[timeIndexMecToCal]
    quatECI2GSEToCal  = quatECI2GSEArr[timeIndexMecToCal]

    quatECI2DBCSToCalQ = pyquaternion.Quaternion(quatECI2DBCSToCal)
    quatECI2GSEToCalQ  = pyquaternion.Quaternion(quatECI2GSEToCal)
    quatECI2DBCSToCalN = quatECI2DBCSToCalQ.normalised
    quatDBCS2ECIToCalQ = quatECI2DBCSToCalN.inverse
    quatDBCS2ECIToCalN = quatDBCS2ECIToCalQ.normalised
    quatECI2GSEToCalN  = quatECI2GSEToCalQ.normalised
    vDBCSVector        = np.hstack((vxKMS, vyKMS, vzKMS))
    vECIVector         = quatDBCS2ECIToCalN.rotate(vDBCSVector)
    vGSEVector         = quatECI2GSEToCalN.rotate(vECIVector)


'''
https://omniweb.gsfc.nasa.gov/ftpbrowser/bow_derivation.html
Converting units, this becomes
        FP = (2*10**-6)*Np*Vp**2 nPa (N in cm**-3, Vp in km/s)
'''
normal = np.array(normal)
vPro   = vGSEVector[0]*normal[0] + vGSEVector[1]*normal[1] + vGSEVector[2]*normal[2]
pdyn   = (2*10**(-6))*densityCM*vPro**2


print('Thermal Pressure of reflected ion beam is', pthnPa, 'nPa')
print('Velocity in GSE is', vGSEVector)
print('Dynamic Pressure of reflected ion beam is', pdyn,   'nPa')

fileNameFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2017/11/16/mms1_fgm_brst_l2_20171116121223_v5.112.0.cdf')
fgmInfo     = fileNameFgm.cdf_info()
epochFgm    = fileNameFgm.varget('Epoch')
bGSEVec     = fileNameFgm.varget('mms1_fgm_b_gse_brst_l2')

epochToCalculateType2000 = cdflib.cdfepoch.compute_tt2000(timeInput)
timeFgmIndexTuple        = np.where(epochFgm > epochToCalculateType2000)
timeFgmIndexArr          = timeFgmIndexTuple[0]
timeFgmIndexToCal        = timeFgmIndexArr[0]

btotalTime = bGSEVec[timeFgmIndexToCal, 3]

pb    = (btotalTime*10**(-9))**2/(2*mu0)
pbnPa = pb*10**9###magnetic pressure in nPa

print('Magnetic Pressure is', pbnPa, 'nPa')

# pyEndTime = time.time()
# print('==================================================================================================================')
# print('==================================================================================================================')
# print('Took %s seconds to run.' % (pyEndTime - pyStartTime))


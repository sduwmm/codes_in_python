import cdflib
import numpy as np
import math
import matplotlib.pyplot as plt

def MomentWMM(Time, FileDist, Ener, Phi, Theta):
    ###constants==================
    mp   = 1.67262158 * 10 ** (-27)
    eV2J = 1.602176462 * 10 ** (-19)
    kb   = 1.3806503 * 10 ** (-23)
    mu0  = 4 * math.pi * 10 ** (-7)
    ###==================constants
    ###init====================
    EnerInput  = np.array(Ener)
    PhiInput   = np.array(Phi)
    ThetaInput = np.array(Theta)
    ###========================init
    ###get data==========================
    infoDist       = FileDist.cdf_info()
    epochDist      = FileDist.varget('epoch')
    phiRaw         = FileDist.varget('mms1_dis_phi_brst')
    phiDeltaRaw    = FileDist.varget('mms1_dis_phi_delta_brst')
    thetaRaw       = FileDist.varget('mms1_dis_theta_brst')
    thetaDeltaRaw  = FileDist.varget('mms1_dis_theta_delta_brst')
    distRaw        = FileDist.varget('mms1_dis_dist_brst')
    energyRaw      = FileDist.varget('mms1_dis_energy_brst')
    energyDeltaRaw = FileDist.varget('mms1_dis_energy_delta_brst')


    epochToCalculateType2000 = cdflib.cdfepoch.compute_tt2000(Time)
    timeDistIndexTuple       = np.where(epochDist > epochToCalculateType2000)
    timeDistIndexArr         = timeDistIndexTuple[0]
    timeDistIndexToCal       = timeDistIndexArr[0]

    ###convert data to new data in SI units==================================
    phi        = phiRaw * math.pi / 180  ###angle in radian
    theta      = thetaRaw * math.pi / 180  ###angle in radian
    deltaPhi   = phiDeltaRaw * math.pi / 180 * 2  ###angle in radian
    deltaTheta = thetaDeltaRaw * math.pi / 180 * 2  ###angle in radian
    dist       = distRaw * 10 ** 12  ###dist in s3/m6
    energy     = energyRaw * eV2J  ###energy in J
    deltaEner  = energyDeltaRaw * eV2J  ###delta energy in J
    ###==================================convert data to new data in SI units

    phiTT       = phi[timeDistIndexToCal]
    thetaTT     = theta
    deltaPhi    = deltaPhi[0]
    deltaTheta  = deltaTheta[0]
    distTT      = dist[timeDistIndexToCal]
    energyTT    = energy[timeDistIndexToCal]
    deltaEnerTT = deltaEner[timeDistIndexToCal]

    LoopEner  = EnerInput
    LoopPhi   = PhiInput
    LoopTheta = ThetaInput

    denTT = 0
    nByVxTT = 0
    nByVyTT = 0
    nByVzTT = 0
    pxxItem1 = 0
    pyyItem1 = 0
    pzzItem1 = 0

    for ii in LoopEner:
        energyII    = energyTT[ii]
        deltaEnerII = deltaEnerTT[ii]
        velII       = math.sqrt(2 * energyII / mp)
        velIISqrt   = velII ** 2
        velIICube   = velII ** 3
        velIIQuar   = velII ** 4
        enerIILeft  = energyII - deltaEnerII
        enerIIRight = energyII + deltaEnerII
        velIILeft   = math.sqrt(2 * enerIILeft / mp)
        velIIRight  = math.sqrt(2 * enerIIRight / mp)
        velDeltaII  = velIIRight - velIILeft
        for jj in LoopTheta:
            thetaJJ = thetaTT[jj]
            for kk in LoopPhi:
                phiKK   = phiTT[kk]
                denTT   += deltaTheta * deltaPhi * velIISqrt * velDeltaII * math.sin(thetaJJ) * distTT[ii, jj, kk]
                nByVxTT += (-1) * deltaTheta * deltaPhi * velIICube * velDeltaII * (math.sin(thetaJJ)) ** 2 * math.cos(
                    phiKK) * distTT[ii, jj, kk]
                nByVyTT += (-1) * deltaTheta * deltaPhi * velIICube * velDeltaII * (math.sin(thetaJJ)) ** 2 * math.sin(
                    phiKK) * distTT[ii, jj, kk]
                nByVzTT += (-1) * deltaTheta * deltaPhi * velIICube * velDeltaII * math.sin(thetaJJ) * math.cos(
                    thetaJJ) * distTT[ii, jj, kk]
                pxxItem1 += mp * deltaTheta * deltaPhi * velIIQuar * velDeltaII * (math.sin(thetaJJ)) ** 3 * (
                    math.cos(phiKK)) ** 2 * distTT[ii, jj, kk]
                pyyItem1 += mp * deltaTheta * deltaPhi * velIIQuar * velDeltaII * (math.sin(thetaJJ)) ** 3 * (
                    math.sin(phiKK)) ** 2 * distTT[ii, jj, kk]
                pzzItem1 += mp * deltaTheta * deltaPhi * velIIQuar * velDeltaII * math.sin(thetaJJ) * (
                    math.cos(thetaJJ)) ** 2 * distTT[ii, jj, kk]
    #
    vxTT = nByVxTT / denTT
    vyTT = nByVyTT / denTT
    vzTT = nByVzTT / denTT
    vtotalTT = math.sqrt(vxTT ** 2 + vyTT ** 2 + vzTT ** 2)
    #
    pxxTT = pxxItem1 - mp * denTT * vxTT ** 2
    pyyTT = pyyItem1 - mp * denTT * vyTT ** 2
    pzzTT = pzzItem1 - mp * denTT * vzTT ** 2
    pTT = (pxxTT + pyyTT + pzzTT) / 3
    #
    vxTT = nByVxTT / denTT
    vyTT = nByVyTT / denTT
    vzTT = nByVzTT / denTT
    vtotalTT = math.sqrt(vxTT ** 2 + vyTT ** 2 + vzTT ** 2)
    #
    pxxTT = pxxItem1 - mp * denTT * vxTT ** 2
    pyyTT = pyyItem1 - mp * denTT * vyTT ** 2
    pzzTT = pzzItem1 - mp * denTT * vzTT ** 2
    pTT = (pxxTT + pyyTT + pzzTT) / 3
    #
    densityM = denTT  ###number density in m**(-3)
    densityCM = denTT * 10 ** (-6)  ###number density in cm*(-3)
    vxMS = vxTT  ###velocity in m/s
    vyMS = vyTT  ###velocity in m/s
    vzMS = vzTT  ###velocity in m/s
    vtotalMS = vtotalTT  ###velocity in m/s
    vxKMS = vxTT * 10 ** (-3)  ###velocity in km/s
    vyKMS = vyTT * 10 ** (-3)  ###velocity in km/s
    vzKMS = vzTT * 10 ** (-3)  ###velocity in km/s
    vtotalKMS = vtotalTT * 10 ** (-3)  ###velocity in km/s
    pthPa = pTT  ###p in Pa
    pthnPa = pthPa * 10 ** 9  ###p in nPa

    tempK = pthPa / (densityM * kb)
    tempeV = (pthPa / densityM) / eV2J

    return(densityCM, vxKMS, vyKMS, vzKMS, tempeV)

###test
TimeInput = [2019, 2, 23, 5, 22, 3, 000]
FileDistInput = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-dist/2019/02/23/mms1_fpi_brst_l2_dis-dist_20190223052123_v3.3.0.cdf')
#EnerInput = [17,26]
#PhiInput = [0,1, 29,  30, 31]
#ThetaInput = [3,11]
EnerInput   = np.arange(32)
PhiInput   = np.arange(32)
ThetaInput = np.arange(16)
moment = MomentWMM(TimeInput,FileDistInput,EnerInput,PhiInput,ThetaInput)






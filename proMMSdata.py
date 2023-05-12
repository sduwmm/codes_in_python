# -*- coding: utf-8 -*-
'''
Created on Tue May 19 16:48:34 2020

@author: Mengmeng Wang
'''

'''
'''

'''
!!!Note
'''

import cdflib
import  numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import time

pyStartTime = time.time()
figurePath = 'C:\\Users\\BRIAR\\Desktop\\OneDrive\\work\\code&figure\\'
if not os.path.exists(figurePath):
    os.makedirs(figurePath)
fileFgm  = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2017/12/01/mms1_fgm_brst_l2_20171201143933_v5.114.0.cdf')
infoFgm  = fileFgm.cdf_info()
epochFgm = fileFgm.varget('Epoch')
bGSEVec  = fileFgm.varget('mms1_fgm_b_gse_brst_l2')
fileDisMom  = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/dis-moms/2017/12/01/mms1_fpi_brst_l2_dis-moms_20171201143933_v3.3.0.cdf')
infoDisMom  = fileDisMom.cdf_info()
epochDisMom = fileDisMom.varget('Epoch')
ionDenVec   = fileDisMom.varget('mms1_dis_numberdensity_brst')
ionTemp0Vec  = fileDisMom.varget('mms1_dis_temppara_brst')
ionTemp1Vec  = fileDisMom.varget('mms1_dis_tempperp_brst')
ionVelVec    = fileDisMom.varget('mms1_dis_bulkv_gse_brst')
fileDesMom  = cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/des-moms/2017/12/01/mms1_fpi_brst_l2_des-moms_20171201143933_v3.3.0.cdf')
infoDesMom  = fileDesMom.cdf_info()
epochDesMom = fileDesMom.varget('Epoch')
eleDenVec   = fileDesMom.varget('mms1_des_numberdensity_brst')
eleTemp0Vec = fileDesMom.varget('mms1_des_temppara_brst')
eleTemp1Vec = fileDesMom.varget('mms1_des_tempperp_brst')


ionTempVec = 1/3*ionTemp0Vec + 2/3*ionTemp1Vec
eleTempVec = 1/3*eleTemp0Vec + 2/3*eleTemp1Vec
bxGSEVec = bGSEVec[:, 0]
byGSEVec = bGSEVec[:, 1]
bzGSEVec = bGSEVec[:, 2]
btGSEVec = bGSEVec[:, 3]

trDown0 = [2017, 12, 1, 14, 40, 56, 500]
trDown1 = [2017, 12, 1, 14, 40, 58, 600]

trUp0   = [2017, 12, 1, 14, 41,  0, 200]
trUp1   = [2017, 12, 1, 14, 41,  5, 000]

#
epoch0DownType2000 = cdflib.cdfepoch.compute_tt2000(trDown0)
epoch1DownType2000 = cdflib.cdfepoch.compute_tt2000(trDown1)
trDownTuple0       = np.where(epochFgm > epoch0DownType2000)
trDownTuple1       = np.where(epochFgm < epoch1DownType2000)
trArrDown0         = trDownTuple0[0]
trArrDown1         = trDownTuple1[0]
trDownArrToCal     = np.intersect1d(trArrDown0, trArrDown1)
#
epochDown          = epochFgm[trDownArrToCal]
bxGSEDown          = bxGSEVec[trDownArrToCal]
byGSEDown          = byGSEVec[trDownArrToCal]
bzGSEDown          = bzGSEVec[trDownArrToCal]
btGSEDown          = btGSEVec[trDownArrToCal]
#
#
epoch0UpType2000 = cdflib.cdfepoch.compute_tt2000(trUp0)
epoch1UpType2000 = cdflib.cdfepoch.compute_tt2000(trUp1)
trUpTuple0       = np.where(epochFgm > epoch0UpType2000)
trUpTuple1       = np.where(epochFgm < epoch1UpType2000)
trArrUp0         = trUpTuple0[0]
trArrUp1         = trUpTuple1[0]
trUpArrToCal     = np.intersect1d(trArrUp0, trArrUp1)
#
epochUp          = epochFgm[trUpArrToCal]
bxGSEUp          = bxGSEVec[trUpArrToCal]
byGSEUp          = byGSEVec[trUpArrToCal]
bzGSEUp          = bzGSEVec[trUpArrToCal]
btGSEUp          = btGSEVec[trUpArrToCal]
#
bxGSEDownMean = np.mean(bxGSEDown)
byGSEDownMean = np.mean(byGSEDown)
bzGSEDownMean = np.mean(bzGSEDown)
bxGSEUpMean   = np.mean(bxGSEUp)
byGSEUpMean   = np.mean(byGSEUp)
bzGSEUpMean   = np.mean(bzGSEUp)

bVectorDown = [bxGSEDownMean, byGSEDownMean, bzGSEDownMean]
bVectorUp   = [bxGSEUpMean, byGSEUpMean, bzGSEUpMean]
bVectorDown = np.array(bVectorDown)
bVectorUp   = np.array(bVectorUp)

numerator    = np.cross(np.cross(bVectorUp, bVectorDown), (bVectorUp - bVectorDown))
denominator0 = np.cross(np.cross(bVectorUp, bVectorDown),(bVectorUp - bVectorDown))
denominator  = math.sqrt(denominator0[0]**2 + denominator0[1]**2 + denominator0[2]**2)
normal = numerator/denominator

print(normal)



fig = plt.figure(figsize = (7,12))
fig.subplots_adjust(top = 0.94)
fig.subplots_adjust(bottom = 0.09)
fig.subplots_adjust(left = 0.13)
fig.subplots_adjust(right = 0.92)
gs = GridSpec(6, 1, figure = fig)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])
ax3 = fig.add_subplot(gs[3])
ax4 = fig.add_subplot(gs[4])
ax5 = fig.add_subplot(gs[5])

ax0.plot(epochDown, bxGSEDown)
ax0.plot(epochDown, btGSEDown, color = 'black')
ax1.plot(epochUp, btGSEUp, color = 'black')
'''
ax0.plot(epochFgm, bGSEVec)
ax1.plot(epochDisMom, ionDenVec)
ax2.plot(epochDisMom, ionTempVec)
ax3.plot(epochDisMom, ionVelVec)
ax4.plot(epochDesMom, eleDenVec)
ax5.plot(epochDesMom, eleTempVec)
'''
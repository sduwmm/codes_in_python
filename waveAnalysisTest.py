
import cdflib
import  numpy as np
import math
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import time

pyStartTime = time.time()
figurePath = 'C:\\Users\\BRIAR\\Desktop\\OneDrive\\work\\code&figure\\'
if not os.path.exists(figurePath):
    os.makedirs(figurePath)
fileFgm    = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2017/12/01/mms1_fgm_brst_l2_20171201143933_v5.114.0.cdf')
infoFgm    = fileFgm.cdf_info()
epochFgm   = fileFgm.varget('Epoch')
bGSEFgmVec = fileFgm.varget('mms1_fgm_b_gse_brst_l2')
bGSEFgm    = bGSEFgmVec[:,0:3]
fileEdp    = cdflib.CDF('D:/data/mms/mms1/edp/brst/l2/dce/2017/12/01/mms1_edp_brst_l2_dce_20171201143933_v3.0.0.cdf')
infoEdp    = fileEdp.cdf_info()
epochEdp   = fileEdp.varget('mms1_edp_epoch_brst_l2')
dceEdp     = fileEdp.varget('mms1_edp_dce_gse_brst_l2')



#
trWaves0 = [2017, 12, 1, 14, 41, 0, 200]
trWaves1 = [2017, 12, 1, 14, 41, 5, 000]
#
epoch0WavesType2000 = cdflib.cdfepoch.compute_tt2000(trWaves0)
epoch1WavesType2000 = cdflib.cdfepoch.compute_tt2000(trWaves1)
#
trFgmTuple0 = np.where(epochFgm > epoch0WavesType2000)
trFgmTuple1 = np.where(epochFgm < epoch1WavesType2000)
trFgmArr0   = trFgmTuple0[0]
trFgmArr1   = trFgmTuple1[0]
trFgmToCal  = np.intersect1d(trFgmArr0, trFgmArr1)
#
epochFgmToCal = epochFgm[trFgmToCal]
#
trEdpTuple0 = np.where(epochEdp > epochFgmToCal[0])
trEdpTuple1 = np.where(epochEdp < epochFgmToCal[-1])
trEdpArr0   = trEdpTuple0[0]
trEdpArr1   = trEdpTuple1[0]
trEdpToCal  = np.intersect1d(trEdpArr0, trEdpArr1)
#
epochEdpToCal = epochEdp[trEdpToCal]
bGSEToCal     = bGSEFgm[trFgmToCal, :]
eDceToCal     = dceEdp[trEdpToCal, :]
#
fx = interpolate.interp1d(epochFgmToCal, bGSEToCal[:,0])
fy = interpolate.interp1d(epochFgmToCal, bGSEToCal[:,1])
fz = interpolate.interp1d(epochFgmToCal, bGSEToCal[:,2])
bxGSEN = fx(epochEdpToCal)
byGSEN = fy(epochEdpToCal)
bzGSEN = fz(epochEdpToCal)
exGSE = eDceToCal[:,0]
eyGSE = eDceToCal[:,1]
ezGSE = eDceToCal[:,2]
#
##transpose
fig = plt.figure(figsize = (7,12))
#fig.subplots_adjust(top = 0.94)
gs = GridSpec(4,1, figure = fig)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])
ax3 = fig.add_subplot(gs[3])
ax0.plot(epochFgmToCal, bGSEToCal[:,0], color = 'blue')
ax0.plot(epochEdpToCal, bxGSEN, color = 'red')
ax1.plot(epochFgmToCal, bGSEToCal[:,1], color = 'blue')
ax1.plot(epochEdpToCal, byGSEN, color = 'red')
ax2.plot(epochFgmToCal, bGSEToCal[:,2], color = 'blue')
ax2.plot(epochEdpToCal, bzGSEN, color = 'red')

'''
epochEdpToCal
exGSE
eyGSE
ezGSE
bxGSE
byGSE
bzGSE
'''
'''
def pf2d0(epochToCal, ex, ey, ez, bx, by, bz):
    noPoints   = len(bx)
    stepLength = 128*2 ###
    nopFFT      = stepLength*2###

    noSteps = (noPoints - nopFFT)/stepLength
    noSteps = int(noSteps)
    nopFFT  = int(nopFFT)
    power   = np.zeros(shape = (noSteps, nopFFT/2))


    return(power)
'''
noPoints = len(bxGSEN)
stepLength = 128*2###
nopFFT     = stepLength*2###

noSteps = (noPoints - nopFFT)/stepLength
noSteps = int(noSteps)
halfNopFFT  = int(nopFFT/2)
power   = np.zeros(shape = (noSteps, halfNopFFT))
bDataT  = np.zeros(shape = (noSteps, 3))
eDataT  = np.zeros(shape = (noSteps, 3))
#eSpecX  = np.zeros(shape = (noSteps, 3)) + np.zeros(shape = (noSteps, 3))j
'''
epochEdpToCal
exGSE
eyGSE
ezGSE
bxGSE
byGSE
bzGSE
'''
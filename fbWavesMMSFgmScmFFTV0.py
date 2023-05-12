
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
bGSEFgm    = fileFgm.varget('mms1_fgm_b_gse_brst_l2')
#bGSEFgmVec = fileFgm.varget('mms1_fgm_b_gse_brst_l2')
#bGSEFgm    = bGSEFgmVec[:,0:3]
fileScm    = cdflib.CDF('D:/data/mms/mms1/scm/brst/l2/scb/2017/12/01/mms1_scm_brst_l2_scb_20171201143933_v2.2.0.cdf')
infoScm    = fileScm.cdf_info()
epochScm   = fileScm.varget('Epoch')
acbGSEScm  = fileScm.varget('mms1_scm_acb_gse_scb_brst_l2')

#
#
#trWaves0 = [2017, 12, 1, 14, 41, 0, 200] #200
#trWaves1 = [2017, 12, 1, 14, 41, 3, 900] #3900
#trWaves0 = [2017, 12, 1, 14, 41, 4, 000]
#trWaves1 = [2017, 12, 1, 14, 41, 5, 000]
#trWaves0 = [2017, 12, 1, 14, 40, 55, 000]
#trWaves1 = [2017, 12, 1, 14, 41, 0, 000 ]
trWaves0 = [2017, 12, 1, 14, 40, 55, 000]
trWaves1 = [2017, 12, 1, 14, 41, 10, 000]
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
bGSEToCal     = bGSEFgm[trFgmToCal, :]
#
#
trScmTuple0 = np.where(epochScm > epoch0WavesType2000)
trScmTuple1 = np.where(epochScm < epoch1WavesType2000)
trScmArr0   = trScmTuple0[0]
trScmArr1   = trScmTuple1[0]
trScmToCal  = np.intersect1d(trScmArr0, trScmArr1)
#
epochScmToCal  = epochScm[trScmToCal]
acbGSEToCal = acbGSEScm[trScmToCal, :]

###FGM FFT
bxForFFT  = bGSEToCal[:,0]
freSam    = 128
tInterval = 1/128
fftSize   = len(epochFgmToCal)
t         = np.arange(0, 1.0, 1.0/freSam)
bxFFT     = np.fft.rfft(bxForFFT)/fftSize
freqs     = np.linspace(0, int(freSam/2), int(fftSize/2+1))
###FGM FFT

###SCM FFT
bxScmForFFT = acbGSEToCal[:,0]
freScmSam   = 8192
fftScmSize  = len(epochScmToCal)
bxFFTScm    = np.fft.rfft(bxScmForFFT)/fftScmSize
freqsScm    = np.linspace(0, int(freScmSam/2), int(fftScmSize/2+1))




fig = plt.figure(figsize = (7,7))
#fig.subplots_adjust(top = 0.94)
gs = GridSpec(4,1, figure = fig)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])
ax3 = fig.add_subplot(gs[3])
#ax0.plot(epochEdpToCal, 1, color = 'blue')
ax0.plot(epochFgmToCal, bGSEToCal)
ax2.plot(epochScmToCal, acbGSEToCal)

ax1.plot(freqs, np.abs(bxFFT))
ax1.set_xlim([0,30])
#ax1.set_ylim([0, 1])

ax3.plot(freqsScm, np.abs(bxFFTScm))
ax3.set_xlim([0,30])
#ax3.set_ylim([0, 0.5])


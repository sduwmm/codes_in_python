import cdflib
import numpy as np
import math
import pyquaternion
import matplotlib.pyplot as plt
import os


timeStartInput = [2017, 12, 1, 14,39, 33, 000]
timeStopInput  = [2017, 12, 1, 14, 42, 43, 000]


epochStartType2000 = cdflib.cdfepoch.compute_tt2000(timeStartInput)
epochStopType2000  = cdflib.cdfepoch.compute_tt2000(timeStopInput)

fileNameFpi= cdflib.CDF('D:/data/mms/mms1/fpi/brst/l2/des-moms/2017/12/01/mms1_fpi_brst_l2_des-moms_20171201143933_v3.3.0.cdf')
infoFpi     = fileNameFpi.cdf_info()
epochFpi    = fileNameFpi.varget('Epoch')
prestensor  = fileNameFpi.varget('mms1_des_prestensor_gse_brst')

timeFpiTuple1    = np.where(epochFpi > epochStartType2000)
timeFpiTuple2    = np.where(epochFpi < epochStopType2000)
timeFpiArr1      = timeFpiTuple1[0]
timeFpiArr2      = timeFpiTuple2[0]
timeFpiArrToPlot = np.intersect1d(timeFpiArr1, timeFpiArr2)
epochFpiToPlot   = epochFpi[timeFpiArrToPlot]

prestensorToPlot=prestensor[timeFpiArrToPlot]
# -*- coding: utf-8 -*-
'''
Created on Wed Jun 5 18:35:30 2019

@author: Mengmeng Wang
'''

'''
function
cut out data in certain time range from cdf object
'''
def cutOutCertainTimeRangeDataFromCdfObj(time1, time2, epochDataType):
    '''
    time1: timeInput1
    time2: timeInput2
    epochDataType: epochFgm
    '''
    epochStartType2000 = cdflib.cdfepoch.compute_tt2000(time1)
    epochStopType2000  = cdflib.cdfepoch.compute_tt2000(time2)
    timeTuple1         = np.where(epochDataType > epochStartType2000)
    timeTuple2         = np.where(epochDataType < epochStopType2000)
    timeArr1           = timeTuple1[0]
    timeArr2           = timeTuple2[0]
    timeArrToCal       = np.intersect1d(timeArr1, timeArr2)
    return(timeArrToCal)

import cdflib
import numpy as np
import math

'''
Call sequence
timeArrFgmToCal = cutOutCertainTimeRangeDataFromCdfObj(timeInput1, timeInput2, epochFgm)
'''


timeInput1 = [2018, 2, 8, 7, 16,  5,000]
timeInput2 = [2018, 2, 8, 7, 17, 15,000]

fileNameFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2018/02/08/mms1_fgm_brst_l2_20180208071603_v5.125.0.cdf')
infoFgm     = fileNameFgm.cdf_info()
epochFgm    = fileNameFgm.varget('Epoch')
bGSEVec     = fileNameFgm.varget('mms1_fgm_b_gse_brst_l2')

timeArrFgmToCal = cutOutCertainTimeRangeDataFromCdfObj(timeInput1, timeInput2, epochFgm)




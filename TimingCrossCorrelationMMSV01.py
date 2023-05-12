'''
Created on Tue Apr 9 11:09:20 2019
Changed on Tue Mar 23 17:33:39 2021
'''
'''
calculate the zero-, first, second- moments 
'''
import cdflib
import numpy as np
import math
#import matplotlib.pyplot as plt
import os
import time

FigurePathInput = 'D:/OneDrive - mail.sdu.edu.cn/work/code&figure/figure2021/skyMap/'

TimeInput = [2017, 12, 18, 12, 56, 48, 500]
Fgm1 = \
    cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2017/12/18/mms1_fgm_brst_l2_20171218125513_v5.117.0.cdf')
#Mec1 =
Fgm2 = \
    cdflib.CDF('D:/data/mms/mms2/fgm/brst/l2/2017/12/18/mms2_fgm_brst_l2_20171218125513_v5.117.0.cdf')

# -*- coding: utf-8 -*-
'''
Created on Mon Jun 04 21:55:23 2019

@author: Mengmeng Wang
'''

import cdflib
import numpy as np
import math
import os
import re

timeInput1  = [2018, 2, 8, 7, 16,  5,000]
timeInput2  = [2018, 2, 8, 7, 17, 15,000]
probes      = ['1']
data_rate_b = 'brst'
level       = 'l2'

dataPath = 'D:/data/mms/'
#fileEndRegex = re.compile(r'(\d)*')
if timeInput1[1] < 10:


    fileFgm = 'D:/data/mms/mms'+str(probes[0])+'/fgm/'+data_rate_b+'/'+level+'/'\
              +str(timeInput1[0])+'/0'+str(timeInput1[1])+'/0'+str(timeInput1[2])+\
              '/mms'+str(probes[0])+'_fgm_'+data_rate_b+'_'+level+'_'+str(timeInput1[0])+'0'+str(timeInput1[1])+'0'+str(timeInput1[2])+\
              '0'+str(timeInput1[3])+'1603_v5.125.0.cdf'

objFgm = cdflib.CDF(fileFgm)

#'D:/data/mms/mms1/fgm/brst/l2/2018/02/08/mms1_fgm_brst_l2_20180208071603_v5.125.0.cdf'
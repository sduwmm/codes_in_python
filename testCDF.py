# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:11:24 2019

@author: WORK
"""

import cdflib
from cdflib import cdfread
from cdflib import cdfwrite
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone, date, time

#datetime.fromtimestamp(timestamp, timezone.utc)


cdfFile = cdflib.CDF('E:/mms1_fgm_brst_l2_20150921102624_v4.18.0.cdf')
info = cdfFile.cdf_info()
print(info)
x = cdfFile.varget('epoch', startrec = 0, endrec = 200)
y = cdfFile.varget("mms1_fgm_b_gsm_brst_l2", startrec = 0, endrec = 200)
print(x)


plt.plot(x,y)

epoch0 = cdflib.cdfepoch.compute_epoch([2015,9,21,10,26,24,111])
print(epoch0)
epoch1 = cdflib.cdfepoch.compute_tt2000([2015,9,21,10,26,24,111])
print(epoch1)
epoch2 = cdflib.cdfepoch.breakdown_tt2000([496103254018796239])
print(epoch2)


tsEpoch = 1362301382
ts = datetime.fromtimestamp(tsEpoch)
ts0 = ts.strftime('%Y-%m-%d %H:%M:%S')
print(ts)

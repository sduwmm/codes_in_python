# -*- coding: utf-8 -*-
'''
Created on Tue Apr 7 15:46:30 2020

@author: Mengmeng Wang
'''

'''
'''

'''
!!!Note
'''

import cdflib
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.gridspec import GridSpec
from numpy import linalg
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import time


timeInput0 = [2017, 12, 1, 14, 40, 00, 000]
timeInput1 = [2017, 12, 1, 14, 41, 30, 000]


###Load FGM Data====================================================================================
fileFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2017/12/01/mms1_fgm_brst_l2_20171201143933_v5.114.0.cdf')
infoFgm     = fileFgm.cdf_info()
epochFgm    = fileFgm.varget('Epoch')

###Load FGM Data====================================================================================


###Load SCM Data====================================================================================
fileScm = cdflib.CDF('D:/data/mms/mms1/scm/brst/l2/scb/2017/12/01/mms1_scm_brst_l2_scb_20171201141123_v2.2.0.cdf')
infoScm     = fileScm.cdf_info()
epochScm    = fileScm.varget('Epoch')
acbGseScm   = fileScm.varget('mms1_scm_acb_gse_scb_brst_l2')
###Load SCM Data====================================================================================
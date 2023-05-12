# -*- coding: utf-8 -*-
'''
Created on Mon Aug 19 16:55:34 2019
Modified on

@author: Mengmeng Wang
'''

import cdflib
import numpy as np
import math
import pyquaternion
import matplotlib.pyplot as plt
from math import pi
#import os




def bowshock(Bz,Mms,Dp,Beta):
    ###constant
    a1 = 11.1266
    a2 = 0.0010
    a3 = - 0.0005
    a4 =  2.5966
    a5 =  0.8182
    a6 = -0.0170
    a7 = -0.0122
    a8 = 1.3007
    a9 = -0.0049
    a10 = -0.0328
    a11 = 6.047
    a12 = 1.029
    a13 = 0.0231
    a14 = -0.002
    epsilon = a12
    ###constant
    if Bz >= 0:
        r0 = a1 *(1 + a2 *Bz) *(1 + a9 *Beta) *(1 + a4 *((a8 - 1) *Mms **2 + 2)/((a8 + 1) *Mms **2)) *Dp **(-1/a11)
        alpha = a5 *(1 + a13 *Bz) *(1 + a7 *Dp) *(1 + a10 *math.log(1 + Beta)) *(1 + a14 * Mms)
    else:
        r0 = a1 *(1 + a3 *Bz) *(1 + a9 *Beta) *(1 + a4 *((a8 - 1) *Mms **2 + 2)/((a8 + 1) *Mms **2)) *Dp **(-1/a11)
        alpha = a5 *(1 + a6 *Bz) *(1 + a7 *Dp) *(1 + a10 *math.log(1 + Beta)) *(1 + a14 *Mms)
    theta1 = np.linspace(0, 0.5 *pi, 1000)
    len_theta1 = len(theta1)
    r1 = np.zeros(len_theta1)
    x1 = np.zeros(len_theta1)
    y1 = np.zeros(len_theta1)
    num_r1 = np.linspace(0, len_theta1 - 1, len_theta1)
    for i in num_r1:
        i = int(i)
        r1[i] = r1[i] + r0 *((1 + epsilon)/(1 + epsilon *np.cos(theta1[i]))) **alpha
        x1[i] = x1[i] + r1[i] *np.cos(theta1[i])
        y1[i] = y1[i] + r1[i] *np.sin(theta1[i])
    theta2 = np.linspace(1.5 *pi, 2 *pi, 1000)
    len_theta2 = len(theta2)
    r2 = np.zeros(len_theta2)
    x2 = np.zeros(len_theta2)
    y2 = np.zeros(len_theta2)
    num_r2 = np.linspace(0, len_theta2 - 1, len_theta2)
    for j in num_r2:
        j = int(j)
        r2[j] = r2[j] + r0 *((1 + epsilon)/(1 + epsilon *np.cos(theta2[j]))) **alpha
        x2[j] = x2[j] + r2[j] *np.cos(theta2[j])
        y2[j] = y2[j] + r2[j] *np.sin(theta2[j])
    x = np.hstack((x1,x2))
    y = np.hstack((y1,y2))
    return x,y
###############################################################################
###############################################################################
Bzn = -0.35
Dpn = 2.48
Mmsn = 6.96
Betan = 2.08
x,y = bowshock(Bzn, Dpn, Mmsn, Betan)
z = y

##########################################################################################

timeStartInput = [2017, 12, 1, 14, 39, 33, 000]
timeStopInput  = [2017, 12, 1, 14, 42, 23, 000]

epochStartType2000 = cdflib.cdfepoch.compute_tt2000(timeStartInput)
epochStopType2000  = cdflib.cdfepoch.compute_tt2000(timeStopInput)

fileNameFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2017/12/01/mms1_fgm_brst_l2_20171201143933_v5.114.0.cdf')
infoFgm     = fileNameFgm.cdf_info()
epochFgm    = fileNameFgm.varget('Epoch')
bGSEVec     = fileNameFgm.varget('mms1_fgm_b_gse_brst_l2')


bxArr = bGSEVec[:,0]
byArr = bGSEVec[:,1]
bzArr = bGSEVec[:,2]
btArr = bGSEVec[:,3]


aaa = np.where(btArr == 34.032116)

bx111   = bxArr[0:10000]
by111   = byArr[0:10000]
bz111   = bzArr[0:10000]
bx1Mean = np.average(bx111)
by1Mean = np.average(by111)
bz1Mean = np.average(bz111)
#[2.0304515, -0.074070625, -0.705954]

bx222   = bxArr[14192:24192]
by222   = byArr[14192:24192]
bz222   = bzArr[14192:24192]
bx2Mean = np.average(bx222)
by2Mean = np.average(by222)
bz2Mean = np.average(bz222)
#[2.4923408, -0.8401598,  -2.2599144]
#bxNew = np.arange(0, bxArr.size)
#byNew = np.arange(0, byArr.size)
#bzNew = np.arange(0, bzArr.size)
#btNew = np.arange(0, btArr.size)

#for i in bxArr.size:


fig = plt.figure(figsize=(16,8))#,dpi=500)
fig.subplots_adjust(top = 0.94)
fig.subplots_adjust(bottom = 0.07)
fig.subplots_adjust(left = 0.13)
fig.subplots_adjust(right = 0.93)

ax1 = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,2)

#ax1.plot(bxArr[:,1], byArr[:,1])



line10 = ax1.scatter(x,y, s = 0.01, color = 'green')
line11 = ax1.scatter(13.2,4.4, s = 10, color = 'black')
ax1.set_xlim([0,20])
ax1.set_ylim([-20,20])



line20 = ax2.scatter(x,z, s = 0.01, color = 'green')
line21 = ax2.scatter(13.2, 2.4, s = 10, color = 'black')
line22 = ax2.plot([13.2,15], [2.4, 15.26], color = 'black')
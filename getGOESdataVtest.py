# -*- coding: utf-8 -*-
'''
Created on Sun Jun 9 21:21:30 2019

@author: Mengmeng Wang
'''
import netCDF4
import numpy as np

ncObj = netCDF4.Dataset('C:/Users/WORK/Desktop/g14_magneto_512ms_20180107_20180107.nc')
var   = ncObj.variables.keys()
time  = ncObj.variables['HE_1'][:,:]
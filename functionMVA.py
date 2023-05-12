# -*- coding: utf-8 -*-
'''
Created on Mom Jun 04 21:18:30 2019

@author: Mengmeng Wang
'''

'''
Minimum Variance Analysis
Magnetic field
Version 0
'''

'''
function MVA:===========================================================================================================
'''
def mva(bx,by,bz):
    mm       = np.ones((3, 3))
    mm[0][0] = np.mean(bx * bx) - np.mean(bx) * np.mean(bx)
    mm[0][1] = np.mean(bx * by) - np.mean(bx) * np.mean(by)
    mm[0][2] = np.mean(bx * bz) - np.mean(bx) * np.mean(bz)
    mm[1][0] = np.mean(by * bx) - np.mean(by) * np.mean(bx)
    mm[1][1] = np.mean(by * by) - np.mean(by) * np.mean(by)
    mm[1][2] = np.mean(by * bz) - np.mean(by) * np.mean(bz)
    mm[2][0] = np.mean(bz * bx) - np.mean(bz) * np.mean(bx)
    mm[2][1] = np.mean(bz * by) - np.mean(bz) * np.mean(by)
    mm[2][2] = np.mean(bz * bz) - np.mean(bz) * np.mean(bz)
    [D, V]  = np.linalg.eig(mm)
    lambda1 = np.min(D)
    Imin    = np.where(D == lambda1)
    lambda3 = np.max(D)
    Imax    = np.where(D == lambda3)
    kLoop   = np.arange(0,3)
    for k in kLoop:
        if D[k] != lambda1 and D[k] != lambda3:
            lambda2 = D[k]
            Imid    = k
    Imin = Imin[0]
    Imin = Imin[0]
    Imax = Imax[0]
    Imax = Imax[0]
    N    = V[:, Imin]
    M    = V[:, Imid]
    L    = V[:, Imax]
    return(lambda1, lambda2, lambda3, N, M, L, V)

'''
!!!Note:
arrange eigenvalues from small to large
lambda1, lambda2, lambda3
corresponding to eigenvectors
N, M, L
corresponding to variation of magnetic field 
minimum, intermediate, maximum
!!!Note
'''
'''
function MVA:===========================================================================================================
'''

'''
Call sequence:
[lambda1, lambda2, lambda3, N, M, L, arrRot] = mva(bxGSEToCal, byGSEToCal, bzGSEToCal)
'''
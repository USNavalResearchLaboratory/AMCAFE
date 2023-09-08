#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 12:09:00 2020

@author: kteferra
"""

import numpy as np
import os
import importlib
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/kteferra/Documents/research/software/PythonLibrary/')
import crystallographicRoutines as cR
import pandas as pd
import time
import h5py

fildir='/Users/kteferra/Documents/research/projects/LLNL'
os.chdir(fildir)

filbaseT='CA3Dnlv8_612'
os.chdir(fildir)
filbase = filbaseT
fil=filbase+'.h5'
f=h5py.File(fil,'r')

nX=f['Step0/dims'][()]
nXs=nX[[2,1,0]] # need to swap x and z b/c x is written out first
dX=f['Step0/VoxelDX'][()]
dXs=dX[[2,1,0]] # need to swap x and z b/c x is written out first
gID=f['Step0/gID'][()]
cTheta=f['Step0/angleAxis'][()]
f.close()
gun=np.unique(gID)
nGrain=gun.shape[0]
aXang=np.reshape(cTheta,[nGrain,4])

vT=np.prod(nX)*np.prod(dX)
nT = np.prod(nX)
volv=np.zeros(nGrain)
for j in range(0,nT):
    volv[gID[j]-1]+=1
wt=volv/nT

# convert ax ang to Euler
eu=cR.AxisAngleToEulerAngle(aXang[:,1:],aXang[:,0])*180.0/np.pi
# get correct grain id's
# write out orientation data
z1=np.concatenate((eu,np.reshape(wt,[nGrain,1])),axis=1)
filout=filbaseT+'_ODF.csv'
np.savetxt(filout,z1,delimiter=",",header="Bunge:phi1,Psi,phi2,Volume %")


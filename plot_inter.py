#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 22:24:50 2020

@author: ogurcan
"""
import numpy as np
import h5py as h5
import matplotlib.pylab as plt

flname='run-GOY/outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
N=kn.shape[0]
Np=11
Sn=np.zeros((Np,N))
for l in range(Np):
    Sn[l,:]=np.mean((np.abs(u[-20000:,:])**2)**(l/2),0)
fl.close()
ftr=np.r_[6:18]
x=np.log10(kn[ftr])
y=np.log10(Sn[:,ftr])
pl=np.zeros((2,Np))
for l in range(Np):
    pl[:,l]=np.polyfit(x,y[l,:],1)
p=np.r_[0:Np]
plt.plot(p,-pl[0,:],'x-',p,p/3,'--')
flname='run-WS04-dynamic/outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
N=kn.shape[0]
Np=11
Sn=np.zeros((Np,N))
for l in range(Np):
    Sn[l,:]=np.mean((np.abs(u[-20000:,:])**2)**(l/2),0)
fl.close()
ftr=np.r_[6:18]
x=np.log10(kn[ftr])
y=np.log10(Sn[:,ftr])
pl=np.zeros((2,Np))
for l in range(Np):
    pl[:,l]=np.polyfit(x,y[l,:],1)
p=np.r_[0:Np]
plt.plot(p,-pl[0,:],'o-')

#flname='outnw.h5'
flname='run-NW04-dynamic/outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
N=kn.shape[0]
Np=11
Sn=np.zeros((Np,N))
for l in range(Np):
    Sn[l,:]=np.mean((np.abs(u[-20000:,:])**2)**(l/2),0)
fl.close()
ftr=np.r_[6:18]
x=np.log10(kn[ftr])
y=np.log10(Sn[:,ftr])
pl=np.zeros((2,Np))
for l in range(Np):
    pl[:,l]=np.polyfit(x,y[l,:],1)
p=np.r_[0:Np]
plt.plot(p,-pl[0,:],'s-')
plt.legend(['GOY','','WS','NW'])

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 09:41:44 2020

@author: ogurcan
"""

import numpy as np
import h5py as h5
import matplotlib.pylab as plt

fig,axes=plt.subplots(5,1,sharex=True,figsize=(8, 10))
#fig.tight_layout()

flname='run-GOY/outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
axes[0].pcolormesh(t[:15000],np.log2(kn),np.log10(np.abs(u[:15000,:])**2+1e-30).T,vmin=-10,vmax=-1,cmap='seismic',rasterized=True)
fl.close()

flname='run-WS04-static/outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
axes[1].pcolormesh(t[:15000],np.log2(kn),np.log10(np.abs(u[:15000,:])**2+1e-30).T,vmin=-10,vmax=-1,cmap='seismic',rasterized=True)
fl.close()

flname='run-NW04-static/outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
axes[2].pcolormesh(t[:15000],np.log2(kn),np.log10(np.abs(u[:15000,:])**2+1e-30).T,vmin=-10,vmax=-1,cmap='seismic',rasterized=True)
fl.close()

flname='run-WS04-dynamic/outnw.h5'
#flname='outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
axes[3].pcolormesh(t[:15000],np.log2(kn),np.log10(np.abs(u[:15000,:])**2+1e-30).T,vmin=-10,vmax=-1,cmap='seismic',rasterized=True)
fl.close()

flname='run-NW04-dynamic/outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
im=axes[4].pcolormesh(t[:15000],np.log2(kn),np.log10(np.abs(u[:15000,:])**2+1e-30).T,vmin=-10,vmax=-1,cmap='seismic',rasterized=True)
plt.xlabel('t')
fl.close()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im,cax=cbar_ax)

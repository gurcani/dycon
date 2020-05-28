#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 16:45:40 2020

@author: ogurcan
"""
import networkx as nx
import h5py as h5
import matplotlib.pylab as plt
import numpy as np


flname='run-NW04-static/outnw.h5'
nwflname='run-NW04-static/nwfile.pkl'
#flname='run-WS04-static/outnw.h5'
#nwflname='run-WS04-static/nwfile.pkl'

fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]

gr=nx.read_gpickle(nwflname)
knl=nx.bipartite.sets(gr)[0]
strs=nx.bipartite.sets(gr)[1]
N=kn.shape[0]
plt.figure(figsize=(6, 4))
plt.loglog(kn,np.mean(np.abs(u[-2000:,:])**2,0)/kn,'x-',kn,1e-4*kn**(-5/3))
plt.axis([1e-1,1e7,1e-20,1e-1])
fl.close()

flname='run-GOY/outnw.h5'
fl=h5.File(flname,"r")
kn=fl['fields/kn'][:]
u=fl['fields/u']
t=fl['fields/t'][:]
plt.loglog(kn,np.mean(np.abs(u[-2000:,:])**2,0)/kn)


pn={kn[l]: np.array([kn[l],1e-17]) for l in range(N)}
pt=dict()
for l in strs:
    exec('a='+l)
    if((a[2]-a[0])==2):
        pt[l]=np.array([kn[a[0]]*2,1e-18])
    else:
        pt[l]=np.array([kn[a[0]]*2,1e-19])
nx.draw_networkx_nodes(kn,pos=pn,node_size=100)
nx.draw_networkx_nodes(strs,pos=pt,node_shape='<',node_size=20)
pos=dict(pn,**pt)
nx.draw_networkx_edges(gr,pos=pos,width=0.5)
ln={kn[l]: l for l in range(N)}
nx.draw_networkx_labels(gr.subgraph(kn),pos=pn,labels=ln,font_size=8,font_color='w')

plt.gca().tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)

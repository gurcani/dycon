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

#nwflname='run-GOY/nwfile.pkl'
#nwflname='run-WS04-static/nwfile.pkl'
nwflname='run-NW04-static/nwfile.pkl'
gr=nx.read_gpickle(nwflname)
kns=nx.bipartite.sets(gr)[0]
kn=np.sort(np.array([l for l in kns]))
strs=nx.bipartite.sets(gr)[1]
N=kn.shape[0]
pn={kn[l]: np.array([np.cos(2*l*np.pi/N),np.sin(2*l*np.pi/N)]) for l in range(N)}
pt=dict()
for l in strs:
    exec('a='+l)
    if((a[2]-a[0])==2):
        pt[l]=0.6*np.array([np.cos(2*(a[0]+1)*np.pi/N),np.sin(2*(a[0]+1)*np.pi/N)])
    else:
        pt[l]=0.4*np.array([np.cos(2*(a[0]+1)*np.pi/N),np.sin(2*(a[0]+1)*np.pi/N)])
plt.figure(figsize=(6, 6))
nx.draw_networkx_nodes(kn,pos=pn,node_size=300)
nx.draw_networkx_nodes(strs,pos=pt,node_shape='<',node_size=20)
pos=dict(pn,**pt)
nx.draw_networkx_edges(gr,pos=pos,width=0.5)
ln={kn[l]: l for l in range(N)}
nx.draw_networkx_labels(gr.subgraph(kn),pos=pn,labels=ln,font_size=10,font_color='w')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 23:40:16 2020

@author: ogurcan
"""

import networkx as nx
import h5py as h5
import matplotlib.pylab as plt
import numpy as np

nwflname='run-GOY/nwfile.pkl'
gr=nx.read_gpickle(nwflname)
cc=nx.bipartite.degrees(gr,nx.bipartite.sets(gr)[0])[1]
kn=np.sort(np.array([l for l in nx.bipartite.sets(gr)[0]]))
N=kn.shape[0]
ccn=[cc[l] for l in kn]
plt.plot(np.arange(N),ccn,'x-')

nwflname='run-WS04-static/nwfile.pkl'
gr=nx.read_gpickle(nwflname)
cc=nx.bipartite.degrees(gr,nx.bipartite.sets(gr)[0])[1]
kn=np.sort(np.array([l for l in nx.bipartite.sets(gr)[0]]))
N=kn.shape[0]
ccn=[cc[l] for l in kn]
plt.plot(np.arange(N),ccn,'s-')

nwflname='run-NW04-static/nwfile.pkl'
gr=nx.read_gpickle(nwflname)
cc=nx.bipartite.degrees(gr,nx.bipartite.sets(gr)[0])[1]
kn=np.sort(np.array([l for l in nx.bipartite.sets(gr)[0]]))
N=kn.shape[0]
ccn=[cc[l] for l in kn]
plt.plot(np.arange(N),ccn,'o-')

plt.legend(['GOY','WS','NW'])
plt.xlabel('n')
plt.ylabel('Degree')

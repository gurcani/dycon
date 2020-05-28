#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 10:47:21 2019

@author: ogurcan
"""
import numpy as np
import time
import scipy.integrate as spi
import h5py as h5
from scipy.stats import norm
import networkx as nx

flname='outnw.h5'

nu=1e-8
t0=0.0                    # initial time
t1=5000                  # final time

dt=1e-3
dtr=1e-2
dtrw=1e-2
dtout=1e-1

N=24
g=2.0
k0=2**(-4)
n=np.arange(0,N)
kn=k0*g**n
trs=[[l,l+1,l+2] for l in range(N-2)]

wecontinue=False            # do we continue from an existing file.
random_forcing=True
dynamic_network=False
save_network=True

#nwstyle='newman-watts'
nwstyle='watts-strogatz'
p=0.0
pf=0.5
def rewire(otrs):
    rtrs=otrs.copy()
    if(nwstyle=='watts-strogatz'):
        for l in range(len(rtrs)):
            if(np.random.rand()<p):
                n=otrs[l][0]
                if(np.random.rand()<pf):
                    if (n<N-3 and n>3):
                        m=np.random.randint(N-n-3)+2
                        rtrs[l]=[n,n+m,n+m+1]
                else:
                    if (n<N-3 and n>3):
                        m=np.random.randint(n-2)+3
                        rtrs[l]=[n,n-m,n-1]
    elif(nwstyle=='newman-watts'):
        for l in range(len(rtrs)):
            if(np.random.rand()<p):
                n=otrs[l][0]
                if(np.random.rand()<pf):
                    if (n<N-3 and n>3):
                        m=np.random.randint(N-n-3)+2
                        rtrs.append([n,n+m,n+m+1])
                else:
                    if (n<N-3 and n>3):
                        m=np.random.randint(n-2)+3
                        rtrs.append([n,n-m,n-1])
    else:
        print('unknown network type')
    return rtrs    

def nonl(l,lp,lpp):
    lp,lpp=np.sort((lp,lpp))
    if(l<lp):
        m=lp-l
        res=kn[lpp]+kn[lp]
    elif(l>lpp):
        m=lpp-lp
        res=((-1)**(lpp-l+2*N)*kn[lpp]-(-1)**(lp-l+2*N)*kn[lp])
    else:
        m=l-lp
        res=((-1)**(lpp-l+2*N)*kn[lpp]-(-1)**(lp-l+2*N)*kn[lp])
    return res/g**m

def connect(otrs):
    lls=[[] for l in range(N)]
    Ms=[[] for l in range(N)]
    for n,l,lp in otrs:
        lls[n].append([l,lp])
        Ms[n].append(nonl(n,l,lp))
        lls[l].append([n,lp])
        Ms[l].append(nonl(l,n,lp))
        lls[lp].append([n,l])
        Ms[lp].append(nonl(lp,n,l))
    ml=0
    for l in lls:
        ml=max(len(l),ml)
    links=np.zeros((N,ml,2),dtype=int)
    Mks=np.zeros((N,ml),dtype=float)
    for l in range(len(lls)):
        for m in range(len(lls[l])):
            links[l,m,:]=lls[l][m]
            Mks[l,m]=Ms[l][m]
    return links,Mks

#
origtrs=trs.copy()
trs=rewire(origtrs)
links,Mks=connect(trs)

Fn=np.zeros(kn.shape,dtype=complex)
Fn[1:3]=1e-2
Dn=nu*kn**2

Gn=Fn.copy()

def force_update():
    global Fn
    Fn=Gn*(norm.rvs(size=(Gn.shape))+1j*norm.rvs(size=(Gn.shape)))

#def func(t,y):
#    u=y.view(dtype=complex)
#    dudt=Fn-Dn*u
#    for n in range(N):
#        for m in range(len(links[n])):
#            dudt[n]+=1j*Mks[n][m]*np.conj(u[links[n][m][0]])*np.conj(u[links[n][m][1]])
#    return dudt.view(dtype=float)


def func(t,y):
    u=y.view(dtype=complex)
    dudt=1j*np.einsum('ij,ij,ij->i',Mks,np.conj(u[links[:,:,0]]),np.conj(u[links[:,:,1]]))+Fn-Dn*u
    return dudt.view(dtype=float)
#
if(wecontinue==True):
    #setting up the output hdf5 file
    fl=h5.File(flname,"r")
    u=fl["fields/u"][:]
    k=fl["fields/k"][:]
    trs=fl["fields/trs"][:]
    tt=fl["fields/t"][:]
    u0=u[-1,:]
    fl.close()
    t0=tt[-1]
    fl=h5.File(flname,"w")
    grp=fl.create_group("fields")
    grp.create_dataset("k",data=k)
    grp.create_dataset("trs",data=trs)
    i=u.shape[0]
    ures=grp.create_dataset("u",(i,N),maxshape=(None,N),dtype=complex)
    tres=grp.create_dataset("t",(i,),maxshape=(None,),dtype=float)
    ures[:,:,:]=u
    tres[:]=tt
else:
#    u0=np.zeros((N),dtype=complex)
#    u0[1:5]=1e-12*np.exp(1j*2*np.pi*np.random.random(N))
    u0=1e-12*np.exp(1j*2*np.pi*np.random.random(N))
    i=0;
    fl=h5.File(flname,"w")
    grp=fl.create_group("fields")
    grp.create_dataset("kn",data=kn)
    ures=grp.create_dataset("u",(1,N),maxshape=(None,N),dtype=complex)
    tres=grp.create_dataset("t",(1,),maxshape=(None,),dtype=float)

#    km0=4*2**(1/4)
#    u0[:,0]=np.exp(-np.linalg.norm(k-k[int(k.shape[0]/2)],axis=0)**2/4**2)*np.exp(1j*np.pi*np.random.random(N))
#    u0[:,1]=np.exp(-np.linalg.norm(k-k[int(k.shape[0]/2)],axis=0)**2/4**2)*np.exp(1j*np.pi*np.random.random(N))
#    u0=np.sqrt(6*np.sqrt(2/np.pi)*km0**(-5)*np.abs(k)**4*np.exp(-2*(np.abs(k)/km0)**2))*np.exp(1j*np.pi*np.random.random(Nh))
if(save_network):
    gr=nx.Graph()
    strs=[np.str(l) for l in trs]
    gr.add_nodes_from(kn,bipartite=0)
    gr.add_nodes_from(strs,bipartite=1)
    
    for l in range(len(trs)):
        gr.add_edges_from([(kn[trs[l][0]] , strs[l]),
                                     (kn[trs[l][1]] , strs[l]),
                                     (kn[trs[l][2]] , strs[l])])
    nx.write_gpickle(gr,'nwfile.pkl')


r=spi.RK45(func,t0,u0.ravel().view(dtype=float),t1,max_step=dt)
epst=1e-12
ct=time.time()
if(random_forcing==True):
    force_update()
#dtff,dtf,dts,dtss=np.sort((dt,dtr,dtrw,dtout))
toldr=-1.0e12
toldrw=-1.0e12
toldout=-1.0e12

while(r.status=='running'):
    told=r.t
    if(r.t>=toldout+dtout-epst and r.status=='running'):    
        toldout=r.t
        print("t=",r.t);
        u=r.y.view(dtype=complex)
        ures.resize((i+1,N))
        tres.resize((i+1,))
        ures[i,:]=u
        tres[i]=r.t
        fl.flush()
        i=i+1;
        print(time.time()-ct,"seconds elapsed.")
    if(r.t>=toldr+dtr-epst and r.status=='running' and random_forcing==True):
        toldr=r.t
        force_update()
    if(r.t>=toldrw+dtrw-epst and r.status=='running' and dynamic_network==True):
        toldrw=r.t
        trs=rewire(origtrs)
        links,Mks=connect(trs)
    while(r.t<told+dt-epst and r.status=='running'):
        res=r.step()
fl.close()

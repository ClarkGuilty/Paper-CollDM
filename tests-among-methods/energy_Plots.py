#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 00:51:41 2019

@author: hiparco
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib import rc, rcParams
from cycler import cycler

rcParams.update({'figure.autolayout': True})

#plt.rcParams['image.cmap'] = 'Blues'
#plt.rcParams['image.cmap'] = 'YlGnBu'
#plt.rcParams['image.cmap'] = 'plasma'
#plt.rcParams['image.cmap'] = 'YlOrRd'
#plt.rcParams['image.cmap'] = 'plasma'
plt.rcParams['image.cmap'] = 'magma'

plt.rcParams['axes.prop_cycle'] = cycler(color='bgrymyk')

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


rc('text', usetex=True)
rcParams.update({'font.size': 16})
fsize = 16

BULLET = -147
JEANS = -137
GAUSS = -127
dt = 0.4


velUnit = 621 #m/s
estUnit = 35 #kpc
potUnit = 385962691092 #J/kg
acceUnit = 3.5737451e-13 #km/s²


constantes = np.loadtxt("constants.dat", usecols = 1)
TAU = int(constantes[8])

dpII = 700
Nt = 100
x = np.linspace(0, Nt*0.1/2, Nt)


condition = '2048-J05-0'



X = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
V = np.linspace(constantes[2], constantes[3], int(constantes[5]))  
#        
f = plt.figure(figsize= (6,5))
t = np.linspace(0, Nt*dt/2,Nt)

f = plt.figure(figsize= (6,5))
#integerLattice = np.loadtxt("2016Method.dat", delimiter=';')

#


for N in 2**np.arange(5,9):

    integerLattice = np.loadtxt('m'+str(N)+"J05.dat", delimiter=';')

    plt.plot(t, integerLattice[:,1]+integerLattice[:,0],  label = str(N))

plt.title("E(t) Mocz casos donde falla Jeans")
plt.xlabel('Tiempo [T]')
plt.ylabel("Energy")
plt.legend()
plt.savefig("JeansFails.png",dpi=dpII)

plt.close(f)
f = plt.figure(figsize= (6,5))

for N in 2**np.arange(9,14):

    integerLattice = np.loadtxt('m'+str(N)+"J05.dat", delimiter=';')

    plt.plot(t, integerLattice[:,1]+integerLattice[:,0],  label = str(N))

plt.title("E(t) Mocz casos donde se activa Jeans")
plt.xlabel('Tiempo [T]')
plt.ylabel("Energy")
plt.legend()
plt.savefig("JeansActivates.png",dpi=dpII)

plt.close(f)
f = plt.figure(figsize= (6,5))

for N in 2**np.arange(5,14):

    integerLattice = np.loadtxt('m'+str(N)+"J05.dat", delimiter=';')

    plt.plot(t, integerLattice[:,1]+integerLattice[:,0],  label = str(N))

plt.title("E(t) Mocz todas las resoluciones")
plt.xlabel('Tiempo [T]')
plt.ylabel("Energy")
plt.legend()
plt.savefig("Jeans.png",dpi=dpII)













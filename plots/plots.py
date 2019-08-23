# -*- coding: utf-8 -*-
"""


@author: Javier Alejandro Acevedo Barroso
Script de Python para la visualización de la simulación.

"""
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib.ticker as ticker
from matplotlib import rc, rcParams
rcParams.update({'figure.autolayout': True})

#plt.rcParams['image.cmap'] = 'Blues'
#plt.rcParams['image.cmap'] = 'YlGnBu'
#plt.rcParams['image.cmap'] = 'plasma'
plt.rcParams['image.cmap'] = 'YlOrRd'
plt.rcParams['image.cmap'] = 'plasma'
#plt.rcParams['image.cmap'] = 'OrRd'

rc('text', usetex=True)
rcParams.update({'font.size': 16})
fsize = 16

BULLET = -147
JEANS = -137
GAUSS = -127
dt = 0.4

#dat = np.loadtxt("./gauss_nocol/grid0.dat").T
##density = np.loadtxt("density.dat")
#constantes = np.loadtxt("./gauss_nocol/constants.dat", usecols = 1)
#TAU = int(constantes[8])

Nt = 100
x = np.linspace(0, 200*0.1/2, 200)

#u0tau0 = np.loadtxt("u0tau0.dat")
#u0tau1 = np.loadtxt("u1tau0.dat")
#plt.figure()
#plt.plot(x,u0tau0, linestyle = ':', label = 'u = 0')
#plt.plot(x,u0tau1, linestyle = '--', label = 'u = $\\sigma$')
#plt.legend()
#plt.ylabel('$\\rho /\\bar{\\rho} - 1$', fontsize = 18)
#plt.xlabel('Time [T]')
#plt.title("Time invariance with $\\tau$ = 0")
#plt.savefig("Jeans_test.png")

plt.figure()
test = np.loadtxt('JeansMagnitude.dat')
test[0] = test[0]/test[0]
plt.plot(x,test)
plt.xlim(0,4)
plt.ylabel('$A_2(t) / A_2(0)$')
plt.yscale('log')


#inF = np.loadtxt("inF.dat")
#outF = np.loadtxt("outF0.dat")
#outF1 = np.loadtxt("outF1.dat")
#oI = np.loadtxt("oI.dat")
#oR = np.loadtxt("oR.dat")
#acce = np.loadtxt("acce.dat")

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

X = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
V = np.linspace(constantes[2], constantes[3], int(constantes[5]))  
#        
constantes[0:4] = constantes[0:4]
figu = plt.gcf()
#figu.set_size_inches(18.5, 10.5)
#figu.set_dpi(300)
dpII = 700
velUnit = 621 #m/s
estUnit = 35 #kpc
potUnit = 385962691092 #J/kg
acceUnit = 3.5737451e-13 #km/s²



condition = 'jeans'
if(condition == 'gauss'): 
      f = plt.figure(figsize= (7,11))
      to_Plot = [49,199]*3
      folders = ['_nocol', '_tau8972', '_tau500']
      ylabels = ['Collisionless', '$\\tau = 8972$', '$\\tau = 500$']
      
      from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
      ax = ImageGrid(f, 111, nrows_ncols=(3,2), axes_pad = 0.2, cbar_mode = 'single', 
                     cbar_pad=0.1, share_all=True,  )
      
      ax[0].set_title(r"t = 5")
      ax[1].set_title(r"t = 20")
      
      for i,j in zip(to_Plot,range(len(to_Plot))):
            dat = np.loadtxt("./"+condition+folders[j//2]+"/grid{:d}.dat".format(i)).T
            
            #ax = f.add_subplot(subplots_index+j)   
            im = ax[j].contourf(np.flip(dat, axis=1), levels = 5, vmin =0,
              extent=[constantes[0],constantes[1],constantes[2],constantes[3]], aspect='auto') #Es mucho más rápido imshow
            ax[j].set_xlim(constantes[0]/2, constantes[1]/2)
            ax[j].set_ylim(constantes[2]/2, constantes[3]/2)
            ax[j].set_ylabel(ylabels[j//2])
            
      ax.cbar_axes[0].colorbar(im)
      ax.cbar_axes[0].set_aspect(30)
      ax.axes_llc.set_xticks([])
      ax.axes_llc.set_yticks([])
      
      f.savefig(condition+'.png', dpi = dpII)






























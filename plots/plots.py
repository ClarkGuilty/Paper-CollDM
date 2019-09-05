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
#plt.rcParams['image.cmap'] = 'YlOrRd'
#plt.rcParams['image.cmap'] = 'plasma'
plt.rcParams['image.cmap'] = 'magma'

rc('text', usetex=True)
rcParams.update({'font.size': 16})
fsize = 16

BULLET = -147
JEANS = -137
GAUSS = -127
dt = 0.4

#dat = np.loadtxt("./gauss_nocol/grid0.dat").T
##density = np.loadtxt("density.dat")
constantes = np.loadtxt("/home/clarkguilty/Paper-CollDM/plots/constants.dat", usecols = 1)
TAU = int(constantes[8])

dpII = 700
Nt = 100
x = np.linspace(0, Nt*0.1/2, Nt)


#plt.figure()
#plt.plot(x,u0tau0, linestyle = ':', label = 'u = 0')
#plt.plot(x,u0tau1, linestyle = '--', label = 'u = $\\sigma$')
#plt.legend()
#plt.ylabel('$\\rho /\\bar{\\rho} - 1$', fontsize = 18)
#plt.xlabel('Time [T]')
#plt.title("Time invariance with $\\tau$ = 0")
#plt.savefig("Jeans_test.png")

plt.figure(figsize=[6.4, 4.8])
#test = np.loadtxt('/home/clarkguilty/Paper-CollDM/plots/JeansMagnitude.dat')
#test[0] = test[0]/test[0]
#plt.plot(x,test)

u0 = np.loadtxt('/home/clarkguilty/Paper-CollDM/plots/invarianceTau0kkj05/u=0.dat')
usigma = np.loadtxt('/home/clarkguilty/Paper-CollDM/plots/invarianceTau0kkj05/u=sigma.dat')
u2sigma = np.loadtxt('/home/clarkguilty/Paper-CollDM/plots/invarianceTau0kkj05/u=2sigma.dat')
u0[0] = u0[0]/u0[0]
usigma[0] = usigma[0]/usigma[0]
u2sigma[0] = u2sigma[0]/u2sigma[0]

#plt.plot(x,u0, color = 'black')
plt.plot(x,u0, color = 'black', label = '$A_2$ with u = 0')
plt.scatter(x[::3],usigma[::3], marker='8', facecolors='red', edgecolors = 'none',label = '$A_2$ with u = $\sigma$')
plt.scatter(x[::3],u2sigma[::3], marker='.', edgecolors = 'green', facecolors= 'none', s=300,label = '$A_2$ with u = $2\sigma$')
plt.xlim(0,5)
plt.ylabel('$A_2(t) / A_2(0)$')
plt.xlabel('Time [T]')
plt.yscale('log')
plt.legend()

plt.title("Perturabación periódica \n$A_2$ vs t con $k/k_j = 0.5$")

plt.savefig('Jeans2Coef.png', dpi =dpII)

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

velUnit = 621 #m/s
estUnit = 35 #kpc
potUnit = 385962691092 #J/kg
acceUnit = 3.5737451e-13 #km/s²



condition = 'gauss'
if(condition == 'gauss'): 
      #f = plt.figure(figsize= (5,11))
      f = plt.figure(figsize= (4,6))
      to_Plot = [49,199]*3
      folders = ['_nocol', '_tau8972', '_tau500']
      ylabels = ['Collisionless', '$\\tau = 8972$', '$\\tau = 500$']
      
      from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
      ax = ImageGrid(f, 111, nrows_ncols=(2,3), axes_pad = 0.2, cbar_mode = 'single', 
                     cbar_pad=0.1, share_all=True, aspect = True, direction= 'column')
      
#      ax[0].set_title(r"t = 5")
#      ax[4].set_title(r"t = 20")
      ax[0].set_ylabel(r"t = 5")
      ax[1].set_ylabel(r"t = 20")


#      ax[0].set_aspect(0.5)
#      ax[1].set_aspect(0.5)
#      ax[2].set_aspect(0.5)
#      ax[3].set_aspect(0.5)
      
      for i,j in zip(to_Plot,range(len(to_Plot))):
            dat = np.loadtxt("./"+condition+folders[j//2]+"/grid{:d}.dat".format(i)).T
            #ax = f.add_subplot(subplots_index+j)   
            im = ax[j].contourf(np.flip(dat, axis=1), levels = 5, vmin =-1,
              extent=[constantes[0]*5.2,constantes[1]*5.2,constantes[2],constantes[3]], aspect='auto') #Es mucho más rápido imshow
            ax[j].set_xlim(constantes[0], constantes[1])
            ax[j].set_ylim(constantes[2]/2, constantes[3]/2)
#            ax[j].set_ylabel(ylabels[j//2]) Descomentar para plot vertical
            ax[j].set_xlabel(ylabels[j//2])
            
      ax.cbar_axes[0].colorbar(im)
      ax.cbar_axes[0].set_aspect(35)
      ax.axes_llc.set_xticks([])
      ax.axes_llc.set_yticks([])
      
      f.savefig(condition+'.png', dpi = dpII, transparent=True)






























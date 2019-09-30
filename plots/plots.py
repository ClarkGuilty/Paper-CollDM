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






X = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
V = np.linspace(constantes[2], constantes[3], int(constantes[5]))  
#        
constantes[0:4] = constantes[0:4]
figu = plt.gcf()



#condition = 'dimensional_invariance05'
#condition = 'gauss'
#condition = 'gauss'
condition = 'galilean_invariance'
t = np.linspace(0, Nt*dt,Nt)

if(condition == 'dimensional_invariance11' ):
    f = plt.figure(figsize= (6,5))
    t = np.linspace(0, Nt*dt,Nt)
#    mass = np.loadtxt(condition+'/massEvolution.dat')[:,0]
    t1024_2_4= np.loadtxt(condition+'/1024-2-4.dat', delimiter=';')
    t2048_4_4= np.loadtxt(condition+'/2048-4-4.dat', delimiter=';')
    t2048_2_1= np.loadtxt(condition+'/2048-2-1.dat', delimiter=';')
    t4096_8_4= np.loadtxt(condition+'/4096-8-4.dat', delimiter=';')
    
    for arr in [t1024_2_4,t2048_4_4,t2048_2_1,t4096_8_4]:
        arr[0] = arr[0]/arr[0]
    
    plt.plot(t, t1024_2_4, label = r'$N_v = 1024$ $k = 2 k_0$ $\bar{\rho} =4$')
    plt.plot(t, t2048_4_4, label = r'$N_v = 2048$ $k = 4 k_0$ $\bar{\rho} =4$')
    plt.plot(t, t2048_2_1, label = r'$N_v = 2048$ $k = 2 k_0$ $\bar{\rho} =1$')
    plt.plot(t, t4096_8_4, label = r'$N_v = 4096$ $k = 8 k_0$ $\bar{\rho} =4$')
    
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.title("Evolución de la energía condición Gaussiana")
    plt.savefig("sadGraph.png", dpi=dpII)    
    
if(condition == 'dimensional_invariance05' ):
    f = plt.figure(figsize= (6,5))
    
#    mass = np.loadtxt(condition+'/massEvolution.dat')[:,0]
    t1024_2_4= np.loadtxt(condition+'/1024-2-4.dat', delimiter=';')
    t2048_4_4= np.loadtxt(condition+'/2048-4-4.dat', delimiter=';')
    t2048_2_1= np.loadtxt(condition+'/2048-2-1.dat', delimiter=';')
    t4096_8_4= np.loadtxt(condition+'/4096-8-4.dat', delimiter=';')
    
    for arr in [t1024_2_4,t2048_4_4,t2048_2_1,t4096_8_4]:
        arr[0] = arr[0]/arr[0]
    
    plt.plot(t, t1024_2_4, label = r'$N_v = 1024$ $k = 2 k_0$ $\bar{\rho} =4$')
    plt.plot(t, t2048_4_4, label = r'$N_v = 2048$ $k = 4 k_0$ $\bar{\rho} =4$')
    plt.plot(t, t2048_2_1, label = r'$N_v = 2048$ $k = 2 k_0$ $\bar{\rho} =1$')
    plt.plot(t, t4096_8_4, label = r'$N_v = 4096$ $k = 8 k_0$ $\bar{\rho} =4$')
    
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.title("Evolución de la energía condición Gaussiana")
    plt.savefig("sadGraph.png", dpi=dpII)    
    


if(condition == 'conservation' ):
    f = plt.figure(figsize= (6,5))
    t = np.linspace(0, Nt*dt/2,Nt)
#    mass = np.loadtxt(condition+'/massEvolution.dat')[:,0]
    energy= np.loadtxt(condition+'/energyEvolution.dat', delimiter=';')
    U0 = energy[0][1]
    plt.plot(t, energy[:,0]-U0, label = r'K')
    plt.plot(t, energy[:,1]-U0, label = r'U')
    plt.plot(t, energy[:,2]-U0, label = r'K+U')
    plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.title("Energy evolution Jeans instability")
    plt.savefig("EvoJeans.png", dpi=dpII)  
#    plt.close(f)
    
    f = plt.figure(figsize= (6,5))
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.title("Energy evolution Jeans instability")
    plt.savefig("EvoJeans.png", dpi=dpII)  
    plt.plot(t, energy[:,4],)
    plt.plot([t[0], t[-1]], [0, 0], color = 'black', label = r'$E_0$')
  
#    plt.close(f)

    

if(condition == 'galilean_invariance'):
    carpetas = ['invarianceTau0kkj05','invarianceTau0kkj11','invarianceTau500kkj05']
    taus = [r'$\infty$',r'$\infty$', '500']
    kas = ['0.5','1.1', '0.5']
    for carpeta, ttau, ka in zip(carpetas,taus,kas):
          
        fig = plt.figure(figsize=[6.4, 4.8])
        u0 = np.loadtxt(carpeta+'/u=0.dat')
        usigma = np.loadtxt(carpeta+'/u=sigma.dat')
        u2sigma = np.loadtxt(carpeta+'/u=2sigma.dat')
#        u0 = u0/u0[0]
#        usigma = usigma/usigma[0]
#        u2sigma = u2sigma/u2sigma[0]
        
        u0[0] = 1
        usigma[0] = 1
        u2sigma[0] = 1
          #plt.plot(x,u0, color = 'black')
        plt.plot(x,u0, color = 'black', label = '$A_2$ with u = 0')
        plt.scatter(x[::3],usigma[::3], marker='8', facecolors='red', edgecolors = 'none',label = '$A_2$ with u = $\sigma$')
        plt.scatter(x[::3],u2sigma[::3], marker='.', edgecolors = 'green', facecolors= 'none', s=150,label = '$A_2$ with u = $2\sigma$')
        plt.xlim(0,5)
        #plt.ylim(1e-4,10)
        plt.ylabel('$A_2(t) / A_2(0)$')
        plt.xlabel('Time [T]')
        plt.yscale('log')
        plt.legend()
         
#        plt.title("Perturabación periódica $\\tau = $"+ttau+"\n$A_2$ vs t con $k/k_j = $"+ka)
        plt.title("Periodic perturbation $\\tau =$ "+ttau+"\n$A_2$ vs t for $k/k_j = $ "+ka)
          
        plt.savefig(carpeta+'/t'+str(ttau)+'Jeans2Coef.png', dpi =dpII)
        plt.close(fig)


if(condition == 'gauss'): 
      #f = plt.figure(figsize= (5,11))
      f = plt.figure(figsize= (6,8))
      to_Plot = [49,199]*3
      folders = ['_nocol', '_tau8972', '_tau500']
      ylabels = ['Collisionless', '$\\tau =$ $8972$', '$\\tau =$ $500$']
      
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






























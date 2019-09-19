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


condition = '2048-J05-0'



X = np.linspace(constantes[0], constantes[1], int(constantes[4]))  
V = np.linspace(constantes[2], constantes[3], int(constantes[5]))  
#        
f = plt.figure(figsize= (6,5))
t = np.linspace(0, Nt*dt/2,Nt)

  
#    mass = np.loadtxt(condition+'/massEvolution.dat')[:,0]
run_char = condition.split(sep='-')
if(run_char[1][0] == 'J'):
    run_char[1] = 'Jeans, $k/k_j$ = '+run_char[1][1]+'.'+run_char[1][2]
energy= np.loadtxt(condition+'/energyEvolution.dat', delimiter=';')
U0 = energy[0][1]
plt.plot(t, energy[:,0]-U0, label = r'K')
plt.plot(t, energy[:,1]-U0, label = r'U')
plt.plot(t, energy[:,2]-U0, label = r'K+U')
plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
plt.legend()
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("Energy evolution\nN = {}, $f_0$ = {}, $\\tau$ = {}".format(run_char[0], run_char[1], run_char[2]))
plt.savefig(condition+"/energy.png", dpi=dpII)  
#    plt.close(f)

f = plt.figure(figsize= (6,5))
plt.xlabel("Time")
plt.ylabel("$1 - E(t)/E(0)$")
plt.title("Energy 'perturbation'")
plt.plot(t, -energy[:,4],)
plt.plot([t[0], t[-1]], [0, 0], color = 'black', label = r'$E_0$')
plt.savefig(condition+"/pertur.png", dpi=dpII)  




















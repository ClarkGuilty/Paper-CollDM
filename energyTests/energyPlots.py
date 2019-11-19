import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sc
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

sep = '-'
#for N in ['512','1024','2048','4096']:
#    for inicial in ['J05']:
#        for tau in ['0','500','8723']:
#            f = plt.figure(figsize= (6,5))            
#            condition = N+sep+inicial+sep+tau
#            print(condition)
#            ax1 = f.add_subplot(111)
##    mass = np.loadtxt(condition+'/massEvolution.dat')[:,0]
#            run_char = condition.split(sep='-')
#            if(run_char[1][0] == 'J'):
#                run_char[1] = 'Jeans, $k/k_j$ = '+run_char[1][1]+'.'+run_char[1][2]
#            energy= np.loadtxt(condition+'/energyEvolution.dat', delimiter=';')
#            U0 = energy[0][1]
#            U0 = 0
##            plt.plot(t, energy[:,0]-U0, label = r'K')
##            U = plt.plot(t, energy[:,1]-U0, label = r'U')
#            E = plt.plot(t, energy[:,2]-U0, label = r'K+U')
#            plt.ylabel("Energy")
#            ax2 = f.add_subplot(111,sharex = ax1, frameon=False)
#            ax2.yaxis.tick_right()
#            ax2.yaxis.set_label_position("right")
#            plt.scatter(t, np.abs((energy[:,2][::3]-U0)/E0) - 1, label = "$|E(t)/E(0)| - 1  $", marker = '+', color = 'r')
##            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
#            handles1, labels1 = ax1.get_legend_handles_labels()
#            handles2, labels2 = ax2.get_legend_handles_labels()
#            handles = handles1+handles2
#            labels = labels1+labels2
#            ax2.legend(handles, labels, loc=3)
#            plt.xlabel("Time")
#            plt.ylabel("$|E(t)/E(0)| - 1  $")
#            plt.title("Energy evolution\nN = {}, $f_0$ = {}, $\\tau$ = {}".format(run_char[0], run_char[1], run_char[2]))
##            plt.savefig(condition+"/energy.png", dpi=dpII)  
#            plt.savefig(condition+".png", dpi=dpII)  
#            plt.close(f)
            
            


#for N in ['512','1024','2048','4096']:
#    for inicial in ['J05', 'J11', 'Gauss']:
#        f = plt.figure(figsize= (6,5))  
#        ax1 = f.add_subplot(111)
#        ax2 = f.add_subplot(111,sharex = ax1, frameon=False)
#        ax2.yaxis.tick_right()
#        ax2.yaxis.set_label_position("right")
#        for tau in ['500','8723','0']:
#            condition = N+sep+inicial+sep+tau
#            print(condition)
#            energy= np.loadtxt(condition+'/energyEvolution.dat', delimiter=';')
#            run_char = condition.split(sep='-')
#            if(tau == '0'): tau = r'$\infty$'
#            if(run_char[1][0] == 'J'):
#                run_char[1] = 'Jeans, $k/k_j$ = '+run_char[1][1]+'.'+run_char[1][2]
#            
#            U0 = energy[0][1]
##            U0 = 0
#            E0 = energy[0][2]-U0
#            ax1.plot(t, energy[:,2]-U0, label = r'$\tau =$ '+tau)
#            
#            handles1, labels1 = ax1.get_legend_handles_labels()
#            handles2, labels2 = [], []
##            ax2.scatter(t[::3], np.abs((energy[:,2][::3]-U0)/E0) - 1, label = r'$\tau =$ '+tau, marker = '+',)
#            if(inicial=='J05' or inicial== 'J11'):
#                ax2.scatter(t[::3], energy[:,4][::3], label = r'$\tau =$ '+tau, marker = '+',)
##            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
#            
#                handles2, labels2 = ax2.get_legend_handles_labels()
#                handles = handles2
#                labels = labels1+labels2
#                ax2.legend(handles, labels, loc=3)
#                ax2.set_ylabel("$|E(t)/E(0)| - 1  $")
#            else:
#                
#                ax1.legend(handles1, labels1, )
#            ax1.set_xlabel("Time")
#            
#            
#            if(inicial[0] == 'G'): plt.axis('off')
#            ax1.set_ylabel("Energy")
#        plt.title("Energy evolution\nN = {}, $f_0$ = {}".format(run_char[0], run_char[1]))
#        plt.savefig('n'+N+'-Evst-'+inicial+".png", dpi=dpII)  
#        plt.close(f)
#




#f = plt.figure(figsize= (6,5))
#integerLattice = np.loadtxt("2016Method.dat", delimiter=';')
#plt.plot(t,integerLattice[:,1]+integerLattice[:,0], label = 'Mocz-Succi')
#
#plt.title('Energy conservation among implementations')
#myIntegerLattice = np.loadtxt("myPoisson.dat", delimiter=';')
#mySpectral = np.loadtxt("mySpectral.dat", delimiter=';')
#mySecond = np.loadtxt("mySecond.dat", delimiter=';')
#plt.ylabel("Energy")
#plt.xlabel("Time [T]")
#plt.plot(t, myIntegerLattice[:,1]+myIntegerLattice[:,0],  label = 'Acevedo-Forero')
##plt.plot(t, mySpectral[:,1]+mySpectral[:,0], label = 'Pseudo-Spectral')
##plt.scatter(t, mySecond[:,1]+mySecond[:,0], label = 'Central-differences', marker ='x')
#plt.legend()
#plt.savefig("energyMocz.png",dpi=dpII)


labels_size = 19

t = np.arange(0,200)*0.1
#for tau in ['0','500']:
#    tau0 = tau
#    for inicial in ['J11']:
#        f = plt.figure(figsize= (6,5))  
#        ax1 = f.add_subplot(111)
##        ax2 = f.add_subplot(111,sharex = ax1, frameon=False)
#        #ax2.yaxis.tick_right()
#        #ax2.yaxis.set_label_position("right")
#        for N1,N2,mark1,mark2,color1,color2,sizes1 in zip(['512','2048'],['1024','4096'],[ 'x','+'],['s','.'],['green','red'],['darkorange','black'],[10,100]):
#            condition1 = N1+sep+inicial+sep+tau0
#            condition2 = N2+sep+inicial+sep+tau0
#            print(condition1,condition2)
#            energy1= np.loadtxt(condition1+'/energyEvolution.dat', delimiter=';')
#            energy2= np.loadtxt(condition2+'/energyEvolution.dat', delimiter=';')
#            for condition in [condition1,condition2]:
#                run_char = condition.split(sep='-')
#                if(tau == '0'): tau = r'$\infty$'
#                if(run_char[1][0] == 'J'):
#                    run_char[1] = 'Jeans, $k/k_j$ = '+run_char[1][1]+'.'+run_char[1][2]
#            
#            U01 = energy1[0][1]
#            U02 = energy2[0][1]
##            U0 = 0
#            E01 = energy1[0][2]-U01
#            E02 = energy2[0][2]-U02
#            #ax1.plot(t, energy[:,2]-U0, label = r'$N =$ '+N)
##            facecolor= color2
##            if(mark1 =='s'):
##                facecolor= 'none'
#            spacing = 12
#
#            ax1.scatter(t[::spacing], np.abs((energy1[:,2][::spacing]-U01)/E01) - 1, label = r'$N =$ '+N1, marker = mark2,edgecolor=color2, facecolor=color2,)
#            ax1.scatter(np.roll(t,spacing//2)[::spacing], np.abs((np.roll(energy2[:,2],3)[::spacing]-U02)/E02) - 1, label = r'$N =$ '+N2, marker = mark1,edgecolor=color1,facecolor=color1)
##            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
#            handles1, labels1 = ax1.get_legend_handles_labels()
##            handles2, labels2 = ax2.get_legend_handles_labels()
#            handles = handles1
#            labels = labels1
#            ax1.legend(handles, labels, )
#            ax1.set_xlabel("Time [T]", fontsize=labels_size)
#            #ax2.set_ylabel("$|E(t)/E(0)| - 1  $")
##            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
##            plt.legend()
##            plt.xlabel("Time")
#            ax1.set_ylabel("$E(t)/E(0)$ - 1",fontsize=labels_size+1)
#        plt.title("Energy evolution\n$\\tau$ = {}, $f_0$ = {}".format(tau, run_char[1]))
#        if(tau == '$\infty$'): tau = '0'
#        plt.savefig('n'+tau0+'-Evst-'+inicial+"-tau"+tau+".png", dpi=dpII)  
#        plt.close(f)



#print('2nd part')
#for inicial in ['J11', 'J05']:
#    for N,in zip(['2048']):
#        f = plt.figure(figsize= (6,5))  
#        ax1 = f.add_subplot(111)
#        #ax2.yaxis.tick_right()
#        #ax2.yaxis.set_label_position("right")
#        for tau,style in zip(['500','10000','0'],['--','-.','-']):    
#            condition = N+sep+inicial+sep+tau
#            print(condition)
#            energy= np.loadtxt(condition+'/energyEvolution.dat', delimiter=';')
#            run_char = condition.split(sep='-')
#            if(tau == '0'): tau = r'$\infty$'
#            if(run_char[1][0] == 'J'):
#                run_char[1] = 'Jeans, $k/k_j$ = '+run_char[1][1]+'.'+run_char[1][2]
#            
#            U0 = energy[0][1]
##            U0 = 0
#            E0 = energy[0][2]-U0
#            #ax1.plot(t, energy[:,2]-U0, label = r'$N =$ '+N)
#            
#
#            ax1.plot(np.linspace(0,np.max(t),len(energy[:,2]))[::3], np.abs((energy[:,2][::3]-U0)/E0) - 1, label = r'$\tau =$ '+tau, linestyle = style)
##            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
#            handles1, labels1 = ax1.get_legend_handles_labels()
##            handles2, labels2 = ax2.get_legend_handles_labels()
#            handles = handles1
#            labels = labels1
#            ax1.legend(handles, labels,fontsize=labels_size-1)
#            plt.xlabel("Time [T]",fontsize=labels_size)
#            #ax2.set_ylabel("$|E(t)/E(0)| - 1  $")
##            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
##            plt.legend()
##            plt.xlabel("Time")
#            ax1.set_ylabel("$E(t)/E(0)$ - 1",fontsize=labels_size+1)
#            toTitle = {'J11' : '1.1', 'J05' : '0.5'}
#        if(tau == '$\infty$'): tau = '0'
#        plt.title("Energy evolution of Jeans instability with \n$k/k_j$ = {} and different relaxation time $\\tau$".format(toTitle[inicial] ))
#        plt.savefig('tau-Evst-'+inicial+".png", dpi=dpII)  
#        plt.close(f)



f = plt.figure(figsize= (8,5))
for implementation, marker in zip(['impl4', 'impl2', 'impl3'],['-','-.', '.']): 
    ax1 = f.add_subplot(111)
    #ax2.yaxis.tick_right()
    #ax2.yaxis.set_label_position("right")
    condition = implementation
    print(condition)
    energy= np.loadtxt(condition+'/energyEvolution.dat', delimiter=';')
    run_char = condition.split(sep='-')
  
    U0 = energy[0][1]
#            U0 = 0
    E0 = energy[0][2]-U0
    #ax1.plot(t, energy[:,2]-U0, label = r'$N =$ '+N)
    
    labels = {'impl1' : 'Acevedo \& Forero + Fuka solver', 
              'impl2' : 'Acevedo \& Forero',
              'impl3' : 'Mocz \& Succi + Fuka',
              'impl4' : 'Acevedo \& Forero +reverse order'}
    ax1.plot(t,energy[:,2]-1, label = labels[implementation])
#            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
    handles1, labels1 = ax1.get_legend_handles_labels()
#            handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles1
    labels = labels1
    ax1.legend(handles, labels,fontsize=labels_size-4, loc = (0.5,0.3))
#    ax1.legend(loc='center left', bbox_to_anchor=(0.5, 0.5))
    plt.xlabel("Time [T]",fontsize=labels_size)
    #ax2.set_ylabel("$|E(t)/E(0)| - 1  $")
#            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
#            plt.legend()
#            plt.xlabel("Time")
#    ax1.set_xlim(0,12)
    ax1.set_ylabel("$E(t)/E(0)$ - 1",fontsize=labels_size+1)
    ax1.set_title("conservación de energía entre implementaciones")
plt.savefig('solvers-Evst-'+".png", dpi=dpII)  
#        plt.close(f)





#print('2nd part')
#for inicial in ['J11', 'J05']:
#    for N,in zip(['2048']):
#        f = plt.figure(figsize= (6,5))  
#        ax1 = f.add_subplot(111)
#        #ax2.yaxis.tick_right()
#        #ax2.yaxis.set_label_position("right")
#        for tau,style in zip(['500','10000','0'],['--','-.','-']):    
#            condition = N+sep+inicial+sep+tau
#            print(condition)
#            energy= np.loadtxt(condition+'/energyEvolution.dat', delimiter=';')
#            run_char = condition.split(sep='-')
#            if(tau == '0'): tau = r'$\infty$'
#            if(run_char[1][0] == 'J'):
#                run_char[1] = 'Jeans, $k/k_j$ = '+run_char[1][1]+'.'+run_char[1][2]
#            
#            U0 = energy[0][1]
##            U0 = 0
#            E0 = energy[0][2]-U0
#            #ax1.plot(t, energy[:,2]-U0, label = r'$N =$ '+N)
#            
#
#            ax1.plot(np.linspace(0,np.max(t),len(energy[:,3]))[::3], np.abs((energy[:,3][::3]-U0)/E0) - 1, label = r'$\tau =$ '+tau, linestyle = style)
##            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
#            handles1, labels1 = ax1.get_legend_handles_labels()
##            handles2, labels2 = ax2.get_legend_handles_labels()
#            handles = handles1
#            labels = labels1
#            ax1.legend(handles, labels,fontsize=labels_size-1)
#            plt.xlabel("Time [T]",fontsize=labels_size)
#            #ax2.set_ylabel("$|E(t)/E(0)| - 1  $")
##            plt.plot([t[0], t[-1]], [energy[0,2], energy[0,2]], color = 'black', label = r'$E_0$')
##            plt.legend()
##            plt.xlabel("Time")
#            ax1.set_ylabel("$E(t)/E(0)$ - 1",fontsize=labels_size+1)
#            toTitle = {'J11' : '1.1', 'J05' : '0.5'}
#        if(tau == '$\infty$'): tau = '0'
#        plt.title("Energy evolution of Jeans instability with \n$k/k_j$ = {} and different relaxation time $\\tau$".format(toTitle[inicial] ))
#        plt.savefig('tau-Evst-'+inicial+".png", dpi=dpII) 


















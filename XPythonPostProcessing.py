import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('C:/Projects/Perry_run')
from origin_write import MySemiconductor

pandas_data = pd.read_csv('C:/Projects/Perry_run/input_csv.csv').to_numpy()
device = MySemiconductor(pandas_data) #put data into MySemiconductor class

project_name = sys.argv[1]
sweeping_v = int(sys.argv[2]) #name of sweeping parameter


def pandas_import_func(filename, Z_Y):
    data = pd.read_csv(filename)
    data = data.to_numpy()
    if Z_Y == 'Y':
        d1 = data[:,0] ; d2 = data[:,5]
    elif Z_Y == 'Z':
        d1 = data[:,0] ; d2 = data[:,6]
    nan = ~pd.isna(d1)
    d1 = d1[nan] ; d2 = d2[nan]
    return( d1, d2)


directory = f'C:/Projects/Projects/{project_name}/CSV'
fig, axs = plt.subplots(1,3, figsize = (12 , 4))

for fname in os.listdir(directory):
    if fname[-3:] == 'csv':
        if fname[-5] =='Y':
            T, Y = pandas_import_func(directory + f'/{fname}' , 'Y')

            axs[0].plot(Y*10**6,T , label = fname)
            if sweeping_v !=0:
                try:
                    axs[2].errorbar(float(fname[:-8]),np.max(T) ,yerr = 2 * np.std(T))
                    axs[2].plot(float(fname[:-8]),np.max(T) , 'ro', markersize=3)
                except:
                    pass
            elif sweeping_v ==0:
                axs[2].errorbar(fname[-7] , np.max(T), yerr = 2 * np.std(T) )
                axs[2].plot(float(fname[-7]),np.max(T) , 'ro', markersize=3)
        elif fname[-5] == 'Z':
            T, Y = pandas_import_func(directory + f'/{fname}' , 'Z')
            axs[1].plot(Y*10**6,T , label = fname)

    else:
        pass

if sweeping_v ==0:
    sweep_title = 'Individual ridges'
elif sweeping_v == 1:
    sweep_title = 'n ridges'
elif sweeping_v ==2:
    sweep_title = 'box x'
elif sweeping_v == 3:
    sweep_title = 'box y'
elif sweeping_v == 4:
    sweep_title = 'box z'
elif sweeping_v == 5:
    sweep_title = 'T sink'
elif sweeping_v == 6:
    sweep_title = 'mesh factor'
elif sweeping_v == 7:
    sweep_title = 'z ridge'
elif sweeping_v == 8:
    sweep_title = 'Power (W)'
else:
    sweep_title = 'no title available'

axs[0].set_title(f'Y sim') ; axs[1].set_title(f'Z sim') ; axs[2].set_title(f'{sweep_title}')
axs[0].set_xlabel('Y - axis (um)') ; axs[0].set_ylabel('Temperature (C)') ; axs[1].set_xlabel('Z - axis (um)') ; axs[1].set_ylabel('Temperature (C)')
axs[0].set_ylim(device.T_sink, ) ; axs[1].set_ylim(device.T_sink, ) ;  axs[2].set_ylim(device.T_sink, )
axs[2].set_ylabel('Temperature (C)') ; axs[2].set_xlabel(f'{sweep_title}')
axs[1].set_xlim(0,)

fig.savefig('plot.png')



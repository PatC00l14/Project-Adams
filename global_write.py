from origin_write import * #import all the functions
import sys
import pandas as pd


def global_writeX(arg0, sweeping_V, V):

    device0 = MySemiconductor(pd.read_csv('input_csv.csv').to_numpy())
    
    if sweeping_V == 8: #sweeping the device power
        for i in range(device0.n_layers):
            device0.r_heat_power[i]  = V * 0.001 / 3 #temporary method to sweep power input
        global_write(arg0 , device0)
    elif sweeping_V ==9:
        device0.thermal_resistance = V
        global_write(arg0 , device0)
    else:
        global_write(arg0 , device0)






    

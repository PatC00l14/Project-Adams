from origin_write import * #import all the functions
import sys
import pandas as pd

device0 = MySemiconductor(pd.read_csv('input_csv.csv').to_numpy())

arg0 = sys.argv[1]
sweeping_V = int(sys.argv[2]) #sweeping argument index
V = float(sys.argv[3]) 

if sweeping_V == 8: #sweeping the device power
    device0.r_heat_power[-1]  = V * 0.001 * (1/8) #temporary method to sweep power input
    global_write(arg0 , device0)
elif sweeping_V ==9:
    V = int(V)
    #device0.r_onoff = np.zeros(device0.n_ridges)
    device0.r_onoff[V] == 1
    global_write(arg0 , device0)
else:
    global_write(arg0 , device0)






    

import numpy as np
import sys 
import os
import shutil
import pandas as pd
sys.path.append('C:/Projects/bin')
from origin_write import MySemiconductor


#Get User input for project name
project_name = input("What is the name of your project?")


print('Which variable would you like to sweep?')


truth_val = True

while True:
    try:
        sweeping_V = int(input("'0 - no sweep \n1 - no. of ridges \n 2 - device x \n 3 - device y \n 4 - device z \n 5 - heat sink temp \n 6 - mesh factor \n 7 - ridge height \n '"))
        break
    except:
        print("You did not enter an integer")


if sweeping_V ==0:
    pass_string = 'ingore'
elif sweeping_V == 1:
    pass_string = 'n_ridges'
elif sweeping_V ==2:
    pass_string = 'box_x'
elif sweeping_V == 3:
    pass_string = 'box_y'
elif sweeping_V == 4:
    pass_string = 'box_z'
elif sweeping_V == 5:
    pass_string = 'T_sink'
elif sweeping_V == 6:
    pass_string = 'mesh_factor'
elif sweeping_V == 7:
    pass_string = 'z_ridge'


#Variable to sweep over
v0 = 200 #initial
v1 = 250 #final
delv = 50 #step


#create results ,logs and input data directory
#Perry_1 -directories and copying of input and processing files
os.chdir(f"C:/Projects/Perry_run")
os.system(f'python Perry_1.py {project_name} {pass_string}')



#os.system(f"python ./global_write.py {project_name}")
#Change working directory to C: drive for Salome
os.system(f"python Perry_active.py {project_name} {v0} {v1} {delv} {sweeping_V}")

#Run Salome to read input CSV and generate the .UNV files stored in an Elmer folder
#After simuation the script writes the final part of the case.sif file with the calculated lowest 
#boundary index
os.chdir(f"../../../Projects/Projects/{project_name}")
os.system(f"python ./XPythonPostProcessing.py")
#Section - 1
#Create the necessary directories, and copy the necessary files to those directories

import numpy as np
import sys 
import os
import shutil
project_name = sys.argv[1] #pass project name as the first and only argument to the file


if not os.path.exists(f'C:/Projects/Projects/{project_name}'): #if already made don't create the new directory
    os.mkdir(f'C:/Projects/Projects/{project_name}')
    os.mkdir(f'C:/Projects/Projects/{project_name}/Logs') ; os.mkdir(f'C:/Projects/Projects/{project_name}/Input_data')
    #create directories for Elmer simulations
    os.mkdir(f'C:/ElmerFEM/ElmerFEM/bin/{project_name}') 
    os.mkdir(f'C:/ElmerFEM/ElmerFEM/bin/{project_name}/dummy')
    os.mkdir(f'C:/Projects/Projects/{project_name}/UNV')
    os.mkdir(f'C:/Projects/Projects/{project_name}/VTU')
    os.mkdir(f'C:/Projects/Projects/{project_name}/CSV')

    print(f"Directory {project_name} created - continue?")
else:
    print("Directory already been used - continue anyway?")

os.system("pause")

#make a copy of the input data  and transfer the python processing script over to the results folder
shutil.copy(f'../bin/Input_csv.csv', f'../Projects//{project_name}/Input_data/input_csv.csv')
shutil.copy(f"../bin/XPythonPostProcessing.py" ,  f"../Projects/{project_name}/XPythonPostProcessing.py")
os.system(f'echo > C:/ElmerFEM/ElmerFEM/bin/{project_name}/case.sif')
#Python script to write up the first part of the case.sif file which Elmer reads to complete
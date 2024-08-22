#Section - 1
#Create the necessary directories, and copy the necessary files to those directories

import numpy as np
import sys 
import os
import shutil
project_name = sys.argv[1] #pass project name as the first and only argument to the file
pass_string = sys.argv[2] #index of sweeping parameter

if not os.path.exists(f'C:/Projects/Projects/{project_name}'): #if already made don't create the new directory
    os.mkdir(f'C:/Projects/Projects/{project_name}')
    os.mkdir(f'C:/Projects/Projects/{project_name}/Logs') ; os.mkdir(f'C:/Projects/Projects/{project_name}/Input_data')
    #create directories for Elmer simulations
    os.mkdir(f'C:/ElmerFEM/ElmerFEM/bin/{project_name}') 
    os.mkdir(f'C:/ElmerFEM/ElmerFEM/bin/{project_name}/dummy')
    os.mkdir(f'C:/Projects/Projects/{project_name}/UNV')
    os.mkdir(f'C:/Projects/Projects/{project_name}/VTU')
    os.mkdir(f'C:/Projects/Projects/{project_name}/CSV')

    while True:
        try:
            print('######################################################################')
            continue_YN = str(input(f"Directory {project_name} created - continue Y/N?"))
            print('######################################################################')
            if continue_YN == 'Y' or continue_YN == 'y':
                break
            elif continue_YN == 'N' or continue_YN == 'n':
                sys.exit()
            else:
                print('Invlaid input. Enter Y,y to continue or N,n to exit')
        except:
            print('some error')

else:
    while True:
        try:
            print('######################################################################')
            print("Directory with this name has already been created - continue anyway? Y/N")
            print('######################################################################')
            continue_YN = str(input('Yy/Nn')) #need option to pull out here
            if continue_YN == 'Y' or continue_YN == 'y':
                break
            elif continue_YN == 'N' or continue_YN == 'n':
                print("Too bad I can\'t seem to terminate the script early \n Just close the terminal and start again")
                break
            else:
                print('Invlaid input. Enter Y,y to continue or N,n to exit')

        except:
            print('some error')



os.system("pause")

#make a copy of the input data  and transfer the python processing script over to the results folder
shutil.copy(f'../bin/Input_csv.csv', f'../Projects//{project_name}/Input_data/input_{pass_string}.csv')
shutil.copy(f"../bin/XPythonPostProcessing.py" ,  f"../Projects/{project_name}/XPythonPostProcessing.py")
shutil.copy(f'Perry_Paraview.py' , f'../../Paraview/Paraview/bin/Perry_Paraview.py')
shutil.copy(f'Perry_Salome.py' , '../../SALOME-9.12.0/W64/Python/Perry_Salome.py')
os.system(f'echo > C:/ElmerFEM/ElmerFEM/bin/{project_name}/case.sif')
#Python script to write up the first part of the case.sif file which Elmer reads to complete
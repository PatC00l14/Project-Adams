#Section - 1
#Create the necessary directories, and copy the necessary files to those directories

import numpy as np
import sys 
import os
import shutil
#project_name = sys.argv[1] #pass project name as the first and only argument to the file
#pass_string = sys.argv[2] #index of sweeping parameter

def create_directories(proj_name):
    """"Create the directories to store input data, simulation files, output data and produce processed plots"""
    if not os.path.exists(f'C:/Projects/Projects/{proj_name}'): #if already made don't create the new directory
        os.mkdir(f'C:/Projects/Projects/{proj_name}')
        os.mkdir(f'C:/Projects/Projects/{proj_name}/Logs') ; os.mkdir(f'C:/Projects/Projects/{proj_name}/Input_data')
        #create directories for Elmer simulations
        os.mkdir(f'C:/ElmerFEM/ElmerFEM/bin/{proj_name}') 
        os.mkdir(f'C:/ElmerFEM/ElmerFEM/bin/{proj_name}/dummy')
        os.mkdir(f'C:/Projects/Projects/{proj_name}/UNV')
        os.mkdir(f'C:/Projects/Projects/{proj_name}/VTU')
        os.mkdir(f'C:/Projects/Projects/{proj_name}/CSV')
        while True:
            try:
                print('######################################################################')
                continue_YN = str(input(f"Directory {proj_name} created - continue Y/N?"))
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

def copy_scripts(proj_name, sweep_ind):
    """"copy the necessary scripts to correct directories for Salome / Elmer / Paraview to operate"""
    shutil.copy(f'C:/Projects/Perry_run/Input_csv.csv', f'C:/Projects/Projects//{proj_name}/Input_data/input_{sweep_ind}.csv')
    #shutil.copy(f'C:/Projects/Perry_run/XPythonPostProcessing.py' ,  f'C:/Projects/Projects/{proj_name}/XPythonPostProcessing.py')
    shutil.copy(f'C:/Projects/Perry_run/Perry_Paraview_nemo.py' , f'C:/Paraview/Paraview/bin/Perry_Paraview_nemo.py')
    shutil.copy(f'C:/Projects/Perry_run/Perry_Salome_nemo.py' , 'C:/SALOME-9.13.0/W64/Python/Perry_Salome_nemo.py')
    shutil.copy(f'C:/Projects/Perry_run/origin_write.py' , 'C:/SALOME-9.13.0/W64/Python/origin_write.py')
    os.system(f'echo > C:/ElmerFEM/ElmerFEM/bin/{proj_name}/case.sif')
    return()

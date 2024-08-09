import numpy as np
import sys 
import os
import shutil
import pandas as pd
sys.path.append('C:/Projects/bin')
from origin_write import MySemiconductor


#Get User input for project name
project_name = input("What is the name of your project?")

#create results ,logs and input data directory
#Perry_1 -directories and copying of input and processing files

os.chdir(f"C:/Projects/Perry_run")
os.system(f'python Perry_1.py {project_name}')

#Variable to sweep over
v0 = 150 #initial
v1 = 300 #final
delv = 50 #step

#os.system(f"python ./global_write.py {project_name}")
#Change working directory to C: drive for Salome
os.system(f"python Perry_active.py {project_name} {v0} {v1} {delv}")

#Run Salome to read input CSV and generate the .UNV files stored in an Elmer folder
#After simuation the script writes the final part of the case.sif file with the calculated lowest 
#boundary index
#redirect to paraview, execute paraview script
os.chdir("C:/Paraview/Paraview/bin")
os.system(f"pvbatch.exe ./Perry_Paraview.py {project_name}")
os.chdir(f"../../../Projects/Projects/{project_name}")
os.system(f"python ./XPythonPostProcessing.py")


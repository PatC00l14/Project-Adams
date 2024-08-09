import pandas as pd
import sys
import shutil
import os 
import numpy as np

sys.path.append('C:/Projects/bin')
from origin_write import MySemiconductor

#import device data using MySemiconductor class - to be called from other programs

project_name = sys.argv[1]
v0 = float(sys.argv[2]) ; v1 = float(sys.argv[3]) ;delv = float(sys.argv[4])

variable = np.arange(v0 , v1, delv)

if True:
    for V in variable:
        #Clear and write bulk of case.sif
        os.chdir('C:/Projects/bin')
        os.system(f"python global_write.py {project_name}")

        #get ready to run salome
        os.chdir('../..')
        os.system(f"SALOME-9.12.0\W64\Python\python3.exe SALOME-9.12.0\salome -t Perry_Salome.py args:{project_name},{V}")

        #get ready for elmer
        os.chdir(f"C:/ElmerFEm/ElmerFEM/bin")
        os.system(f'elmergrid 8 2 {project_name}/temp_save.unv -autoclean -out {project_name}/dummy')
        os.system(f'elmersolver > {project_name}/convergence_log.log 2>&1')

        #copying and saving of design, results and logs to relative folders
        shutil.copy(f"{project_name}/temp_save.unv" , f"../../../Projects/Projects/{project_name}/UNV/{V}.unv")
        shutil.copy(f"{project_name}/case_t0001.vtu" , f"../../../Projects/Projects/{project_name}/VTU/{V}.vtu")
        shutil.copy(f"{project_name}/convergence_log.log" , f"../../../Projects/Projects/{project_name}/Logs/{V}.log")
        #Run Paraview to process output data
        os.chdir("C:/Paraview/Paraview/bin")
        os.system(f"pvbatch.exe ./Perry_Paraview.py {project_name} {V}")
        os.chdir("../../..")
else:
    pass



        




import pandas as pd
import sys
import shutil
import os 
import numpy as np
from global_write import global_writeX
from origin_write import MySemiconductor, global_write1, global_write2
#import device data using MySemiconductor class - to be called from other programs

def perry_active(project_name, v0, v1, delv, sweeping_V):

    #project_name = sys.argv[1]
    #v0 = float(sys.argv[2]) ; v1 = float(sys.argv[3]) ;delv = float(sys.argv[4]
    #sweeping_V = int(sys.argv[5]) #getting sweeping parameter index
    
    if delv ==0: #when the simulation is not sweeping a variable
        variable = np.array([0])
    else:
        variable = np.arange(v0 , v1, delv)
    
    if True:
        for V in variable:
            #Clear and write bulk of case.sif
            os.chdir('C:/Projects/Perry_run')
            #global_writeX(project_name, sweeping_V, V)
    
            #get ready to run salome
            os.chdir('../..')
            os.system(f"SALOME-9.13.0\W64\Python\python3.exe SALOME-9.13.0\salome -t Perry_Salome_nemo.py args:{project_name},{V},{sweeping_V}")
            #get ready for elmer
            os.chdir(f"C:/ElmerFEm/ElmerFEM/bin")
            if not os.path.exists(f'C:/ElmerFEm/ElmerFEM/bin/{proj_name}/temp_save.unv'): 
                os.system(f'elmergrid 8 2 {project_name}/temp_save.unv -autoclean -relh 0.1 -out {project_name}/dummy')
                os.system(f'elmersolver > {project_name}/convergence_log.log 2>&1')
            else:
                sys.exit()
    
            #copying and saving of design, results and logs to relative folders
            shutil.move(f"{project_name}/temp_save.unv" , f"../../../Projects/Projects/{project_name}/UNV/{V}.unv")
            shutil.copy(f"{project_name}/convergence_log.log" , f"../../../Projects/Projects/{project_name}/Logs/{V}.log")
            while True:
                try:
                    shutil.copy(f"{project_name}/case_t0001.vtu" , f"../../../Projects/Projects/{project_name}/VTU/{V}.vtu")
                    os.chdir("C:/Paraview/Paraview/bin")
                    os.system(f"pvbatch.exe ./Perry_Paraview_nemo.py {project_name} {V} {sweeping_V}")
                    os.chdir("../../..")
                except:
                    print("####System likely diverged, cannot solve using Elmer ####")
                break
            #Run Paraview to process output data



        




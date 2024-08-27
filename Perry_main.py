import numpy as np
import sys 
import os
import shutil
import pandas as pd
sys.path.append('C:/Projects/bin')




def main():
    #Get User input for project name to be used to create directories and results locations
    project_name = input("What is the name of your project?")
    
    print('Which variable would you like to sweep?')

    while True:
        try:
            sweeping_V = int(input("'\n0 - no sweep \n1 - no. of ridges \n 2 - device x \n 3 - device y \n 4 - device z \n 5 - heat sink temp \n 6 - mesh factor \n 7 - ridge height \n 8 - Power Sweep \n 9 - arbitrary sweep \n 10 - Device index on/off \n '"))
            break
        except:
            print("You did not enter an integer")
    
    if sweeping_V != 0:
        while True:
            try:
                v0 = float(input("What is starting value of sweep?"))
                print(v0)
                v1 = float(input("What is finishing value of sweep?"))
                print(v1)
                delv = float(input("What is step value of sweep?"))
                print(delv)

                if v0<v1 and delv < np.abs(v1-v0):
                    break
                elif v1 < v0:
                    print('Starting value must be greater than final value')
                elif delv > np.abs(v1-v0):
                    print('Step value is larger than interval space')
                else:
                    print('I have no idea how on earth this error message could\'ve been reached')
                
            except:
                print("Enter a valid value")
    else:
        v0 , v1 , delv = 0 ,0 ,0 #default values if no sweep is chose
    
    
    if sweeping_V ==0:
        pass_string = 'ignore'
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
    elif sweeping_V == 8:
        pass_string = 'Power'
    elif sweeping_V == 9:
        pass_string = 'Arb_sweep'
    elif sweeping_V == 10:
        pass_string = 'Device index on/off'
    
    #Variable to sweep over
    #v0 = 100 #initial
    #v1 = 400 #final
    #delv = 10 #step
    
    
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
    os.system(f"python ./XPythonPostProcessing.py {sweeping_V}")


    print('########################################################')
    print(f'Project: {project_name} complete')
    print('########################################################')


    return(0)


if __name__ == '__main__':
    main()
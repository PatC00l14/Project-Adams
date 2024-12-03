README // Quick explanation of all the python files:

__init__.py:

This is just to initialize the scripts for execution via terminal

Perry_main.py:
Houses the __main__ function. When starting the software direct to this directory and use the command:
"python.exe Perry_main.py" to get everything started
Takes in user project name and sweep parameters then calls the other main functions in the other scripts.

Perry_1.py:
Function to create the necessary results and simulation directories required for the other software, also creates the results directory and creates a copy of the user input parameters. 

Perry_active.py:
Function to actually execute the necessary software while also passing the necessary user input information, such as sweeping variables and the project name. This just uses os. modules to execute commands on terminal so that we can keep everything on python and avoid tedious .bat scripts.

Origin_write.py:
MySemiconductor - this is a class which gets imported by a few other scripts. An easy way to store all the user input information in a simple instance with easy to call attributes. 

Includes all the functions required to converts the user input data into strings and then to write to the case.sif file which is required for ElmerFEM to execute the Multiphysics solvers. 
Since some of the case.sif needs to be written before Salome, some during and the rest after many of the functions are collected into a batch function which then get executed directly from Perry_Salome_nemo.py

Global_write.py:
Partly redundant but previously was the overall method to write the case.sif file. 

Perry_salome_nemo.py:
The python script required to automate the whole meshing and creation of the geometry in Salome. Once perry_1 get activated this gets copied to the C:/Salome-9.13.0/W64/Python directory such that it can interact with the specific Python environment. Since it is a copy that is executed during the run you can continue to make adjustments or code the original version while simulations are taking place.

Perry_paraview_nemo.py:
The exact same as above, instead in the C:/Paraview/paraview/bin directory and for the paraview software

XPythonPostprocessing.py:
Very simply data analysis script, requires specific file names created from the Perry_paraview_nemo.py script. Automatically saves a plot of the .csv files in the results directory.

Salome_internal.py:
This is just a copy of Perry_salome_nemo.py however with the passed arguments removed. This is for debugging only while in the Salome software. This can be used in Salome by selecting 'Load Script'. Helpful to see visually and inspect meshes when there are issues with the Salome side. 



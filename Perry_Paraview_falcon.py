from paraview.simple import *
import numpy as np
import os 
import sys
import pandas as pd
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


sys.path.append('C:/Projects/bin')
from origin_write import MySemiconductor #import semiconductor class



device = pd.read_csv('C:/Projects/bin/input_csv.csv').to_numpy()
device = MySemiconductor(device)




arg0  = sys.argv[1]  #get paraview executable as file/project name
V  = 5
V = float(sys.argv[2])
sweeping_V = int(sys.argv[3])#getting sweeping parameter which is being varied to adjust measurement locations

if sweeping_V == 1:
    device.n_ridges = V
    pass_string = 'n_ridges'
elif sweeping_V ==2:
    device.device_dim[0] = V
    pass_string = 'box_x'
elif sweeping_V == 3:
    device.device_dim[1] = V
    pass_string = 'box_y'
elif sweeping_V == 4:
    device.device_dim[2] = V
    pass_string = 'box_z'
elif sweeping_V == 5:
    device.T_sink = V
    pass_string = 'T_sink'
elif sweeping_V == 6:
    pass_string = 'mesh_factor'
    device.bdy_mesh = device.bdy_mesh *V
    device.r_mesh = device.r_mesh *V
elif sweeping_V == 7:
    device. z_ridge = V
    pass_string = 'z_ridge'
elif sweeping_V ==0:
    pass


# create a new 'XML Unstructured Grid Reader'


def save_data(f_name):

    solved_object = XMLUnstructuredGridReader(registrationName= f'{f_name[-8:]}' , FileName = [f'{f_name}'])
    # set active source
    SetActiveSource(solved_object)
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    
    # show data in view
    w125mesh004vtuDisplay = Show(solved_object, renderView1, 'UnstructuredGridRepresentation')
    
    # trace defaults for the display properties.
    w125mesh004vtuDisplay.Representation = 'Surface'
    
    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    
    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)
    
    # set scalar coloring
    ColorBy(w125mesh004vtuDisplay, ('POINTS', 'temperature'))
    
    # rescale color and/or opacity maps used to include current data range
    w125mesh004vtuDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    w125mesh004vtuDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get 2D transfer function for 'temperature'
    temperatureTF2D = GetTransferFunction2D('temperature')
    
    # get color transfer function/color map for 'temperature'
    temperatureLUT = GetColorTransferFunction('temperature')
    temperatureLUT.TransferFunction2D = temperatureTF2D
    temperatureLUT.RGBPoints = [24.999999999993303, 0.231373, 0.298039, 0.752941, 106.00978698337201, 0.865003, 0.865003, 0.865003, 187.0195739667507, 0.705882, 0.0156863, 0.14902]
    temperatureLUT.ScalarRangeInitialized = 1.0
    
    # get opacity transfer function/opacity map for 'temperature'
    temperaturePWF = GetOpacityTransferFunction('temperature')
    temperaturePWF.Points = [24.999999999993303, 0.0, 0.5, 0.0, 187.0195739667507, 1.0, 0.5, 0.0]
    temperaturePWF.ScalarRangeInitialized = 1
    
    # Properties modified on w125mesh004vtu
    solved_object.TimeArray = 'None'
    
    # show data in view
    w125mesh004vtuDisplay = Show(solved_object, renderView1, 'UnstructuredGridRepresentation')
    
    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)
    
    # show color bar/color legend
    w125mesh004vtuDisplay.SetScalarBarVisibility(renderView1, True)
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # create a new 'Plot Over Line'
    
    
    # Properties modified on plotOverLine1
    
    LinesY = np.array([])
    LinesZ = np.array([])

    [dx ,dy ,dz] = device.device_dim *10**-6
    dz = 99.5*10**-6

    dx = dx - 50 *10**-6

    for i in range(0, device.n_ridges):
        
        LinesY = np.append(LinesY , PlotOverLine(registrationName=f'PlotOverLine{i}', Input=solved_object))
        LinesY[i].Point1 = [dx + (i * 50*10**-6) , 0 , dz] ; LinesY[i].Point2 = [dx + (i * 50*10**-6) , dy, dz]
        SaveData(f'C:/Projects/Projects/{arg0}/CSV/{V}_{i}_Y.csv', LinesY[i], PointDataArrays=['temperature'])


        if device.ext_sink_mat == 0:
            LinesZ = np.append(LinesZ , PlotOverLine(registrationName=f'PlotOverLine{i}', Input=solved_object))        
            LinesZ[i].Point1 = [dx + (i * 50*10**-6) , 0 , 0] ; LinesZ[i].Point2 = [dx + (i * 50*10**-6) , 0, dz]      
        else:
            LinesZ = np.append(LinesZ , PlotOverLine(registrationName=f'PlotOverLine{i}', Input=solved_object))        
            LinesZ[i].Point1 = [dx + (i * 50*10**-6) , 0 , 0] ; LinesZ[i].Point2 = [dx + (i * 50*10**-6) , 0, dz]  
        
        SaveData(f'C:/Projects/Projects/{arg0}/CSV/{V}_{i}_Z.csv', LinesZ[i], PointDataArrays=['temperature'])


    return(0)

directory = f'C:/ElmerFEM/ElmerFEM/bin/{arg0}'
L_dir = len(directory)

f = os.path.join(directory, "case_t0001.vtu")
save_data(f)
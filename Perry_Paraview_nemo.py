from paraview.simple import *
import numpy as np
import os 
import sys
import pandas as pd
sys.path.append('C:/Projects/Perry_run')
from origin_write import MySemiconductor #import semiconductor class

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

device1 = pd.read_csv('C:/Projects/Perry_run/input_csv.csv').to_numpy()
device1 = MySemiconductor(device1)

arg0  = sys.argv[1]  #get paraview executable as file/project name
V = float(sys.argv[2])
sweeping_V = int(sys.argv[3])#getting sweeping parameter which is being varied to adjust measurement locations

if sweeping_V == 1:
    device1.n_ridges = V
    pass_string = 'n_ridges'
elif sweeping_V ==2:
    device1.device_dim[0] = V
    pass_string = 'box_x'
elif sweeping_V == 3:
    device1.device_dim[1] = V
    pass_string = 'box_y'
elif sweeping_V == 4:
    device1.device_dim[2] = V
    pass_string = 'box_z'
elif sweeping_V == 5:
    device1.T_sink = V
    pass_string = 'T_sink'
elif sweeping_V == 6:
    pass_string = 'mesh_factor'
    device1.bdy_mesh = device1.bdy_mesh *V
    device1.r_mesh = device1.r_mesh *V
elif sweeping_V == 7:
    device1. z_ridge = V
    pass_string = 'z_ridge'
elif sweeping_V == 9:
    pass_string = 'Sink_Height'
elif sweeping_V ==0:
    pass

def save_line_data( X1 , X2, V , i , Y_Z, solved_object, c = -1 ):
    [x1,y1,z1] = X1*10**-6 ; [x2, y2, z2] = X2*10**-6 
    #X1,X2:begining / end coordinates. These should be in ""um""
    #V: sweeping variable value. V = 0.0 if no sweep
    #i: laser index number
    #Y_Z: Y axis or Z axis specific
    if c!=-1:
        line_data = PlotOverLine(registrationName=f'PlotOverLine{i}', Input=solved_object)
        line_data.Point1 = [x1 , y1 , z1] ; line_data.Point2 = [x2 , y2, z2]
        SaveData(f'C:/Projects/Projects/{arg0}/CSV/{V}_{c}_{i}_{Y_Z}.csv', line_data, PointDataArrays=['temperature'])
    else:
        line_data = PlotOverLine(registrationName=f'PlotOverLine{i}', Input=solved_object)
        line_data.Point1 = [x1 , y1 , z1] ; line_data.Point2 = [x2 , y2, z2]
        SaveData(f'C:/Projects/Projects/{arg0}/CSV/{V}_0_{i}_{Y_Z}.csv', line_data, PointDataArrays=['temperature'])
    return()

def save_data(f_name,device):
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
    

    [dx ,dy ,dz] = device.device_dim #chip dimensions in um
    dy = dy - 250
    dz = dz + device.z_ridge
    subm_z = device.ext_sink_dim[2] #z scale of submount / heat sink
    n_c = device.n_chips #number of chip
    n_r = device.n_ridges #number of ridges per chip
    temp_array =np.array([1,10,20])
    if device.cartridge_mat !=0: #chuch and cartridge present
        for c in temp_array:
            for r in range(n_r):
                x_pos = (61 + 8.9*(c-1/2))*1000 + dx*(r + 0.5)
                x_1 = np.array([x_pos, 2*1000     , (20 + 12)*1000 + 50 + subm_z + dz])
                x_2 = np.array([x_pos, 2*1000 + dy, (20 + 12)*1000 + 50 + subm_z + dz])
                save_line_data(x_1, x_2, V, r, 'Y',solved_object, c=c)
                x_1 = np.array([x_pos, 2*1000     , (20 + 12)*1000 + 50 + subm_z])
                x_2 = np.array([x_pos, 2*1000 + dy, (20 + 12)*1000 + 50 + subm_z])
                save_line_data(x_1, x_2, V+0.14, r, 'Y',solved_object, c=c)

                x_1 = np.array([x_pos, 2*1000, (20 + 12)*1000 + subm_z+ 50])
                x_2 = np.array([x_pos, 2*1000, (20 + 12)*1000 + subm_z + 50 + dz])
                save_line_data(x_1, x_2, V, r, 'Z',solved_object, c=c)
    else: #Just submount and a single chip
        for r in range(n_r):
            x_pos = dx*(r + 0.5)
            x_1 = np.array([x_pos,  0, dz])
            x_2 = np.array([x_pos, dy, dz])
            save_line_data(x_1, x_2, V, r, 'Y',solved_object)
            x_1 = np.array([x_pos, 0, 0])
            x_2 = np.array([x_pos, 0, dz])
            save_line_data(x_1, x_2, V, r, 'Z',solved_object)
            if device.thermistor_mat != 0:
                tx, ty, tz = device.thermistor_dim
                x_pos += 2 * dx + tx/2
                x_1 = np.array([x_pos, ty/2, 0 ])
                x_2 = np.array([x_pos, ty/2, tz])
                save_line_data(x_1, x_2, V+0.14, r, 'Z',solved_object)

    return(0)

directory = f'C:/ElmerFEM/ElmerFEM/bin/{arg0}'
L_dir = len(directory)

f = os.path.join(directory, "case_t0001.vtu")
save_data(f, device1)
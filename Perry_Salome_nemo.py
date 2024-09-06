import sys
import salome
import numpy as np
import pandas as pd
import os.path

sys.path.append('C:/Projects/Perry_run')
from origin_write import MySemiconductor

#from Perry_active import device_hold
salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0,r'C:/Users/pmillar/Documents/Internship 2024/Week 4/Salome')

import GEOM
from salome.geom import geomBuilder
import SALOMEDS
import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

def write_boundary_conds(boundary1 , arg0 , T_sink = 25): #heat sink temperature - though is not necessarily that complicated to change the input into an arry
    #here we need the boundary number of the lowest boundary (which was far too much effort to determing from the stupid Salome simulations)
    bound_cond = f'\n\nBoundary Condition 1\n  Target Boundaries(1) = {boundary1}\n  Name = "Heat Sink"\n  Temperature = {T_sink}\nEnd'
    mypath = f'C:/ElmerFEM/ElmerFEM/bin/{arg0}/'
    path = os.path.join(mypath , 'case.sif')

    file = open(path , 'a') ; file.write(bound_cond) ; file.close()
    return()


def create_ridges(sc_data , geompy): #this function creates the combined ridge structure
    
    base_y = sc_data.device_dim[1]
    widths = sc_data.r_widths
    thicknesses = sc_data.r_heights
    lengths = sc_data.r_lengths
    n_ridges = sc_data.n_layers
    ridge_partition = np.array([])

    if len(widths) != n_ridges or len(thicknesses)!= n_ridges or len(lengths) != n_ridges:
        raise Exception('Number of layers in input does not match number of layers input')
    
    for i in range(0 , n_ridges):
        ridge = geompy.MakeBoxDXDYDZ(widths[i] , lengths[i] , thicknesses[i])
        ridge = geompy.MakeTranslation(ridge , -widths[i]/2 , (base_y - lengths[i])/2 , np.sum(thicknesses[:i]))
        geompy.addToStudy(ridge , f'ridge{i}')
        ridge_partition = np.append(ridge_partition , ridge)

    combined_ridge = geompy.MakePartition(ridge_partition.tolist(), [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

    return(combined_ridge)


def create_device ( data, geompy, back_facet_trench = False):
    #this will create the a single chip with a submount. The ridge is created within this function

    ridge = create_ridges(data , geompy)

    partition_array = []

    [base_x , base_y , base_z] = data.device_dim
    [trench_x, trench_y , trench_z] = data.trench_dim
    [ext_sink_x,ext_sink_y,ext_sink_z] = data.ext_sink_dim
    n_active = data.n_ridges
    z_pos = base_z  + data.z_ridge

    trench = geompy.MakeBoxDXDYDZ(trench_x , trench_y , trench_z )
    trench = geompy.MakeTranslation(trench , -trench_x , 0 , 0 )

    base = geompy.MakeBoxDXDYDZ(base_x * (n_active + 1) , base_y , base_z)
    base = geompy.MakeTranslation(base , -base_x / 2 , 0 , 0 )

    for i in range( 0 , n_active): #insert all the active regions into the chip
        x_pos = base_x * ( i + 0.5)
        ridge_cut = geompy.MakeTranslation(ridge , x_pos , 0 , z_pos)
        base = geompy.MakeCut(base , ridge_cut , True)
        partition_array = np.append(partition_array, ridge_cut)

    
    for i in range ( 1 , n_active): #create trenches halfway between each active region
        x_pos = base_x * i
        trench_cut = geompy.MakeTranslation(trench , x_pos , 0 , base_z - trench_z)
        base = geompy.MakeCut(base , trench_cut , True)

    if back_facet_trench: #ability to create a back facet trench if there is a back facet monitor present
        back_trench = geompy.MakeBoxDXDYDZ(base_x , trench_x , trench_z)
        back_facet_trench = geompy.MakeTranslation(back_trench, 0 , base_y-330 , base_z - trench_z )
        base = geompy.MakeCut(base , back_facet_trench , True)

    partition_array = np.append(partition_array  ,base) 

    #add submount or external_heatsink
    ext_sink = geompy.MakeBoxDXDYDZ(ext_sink_x,ext_sink_y,ext_sink_z) 
    ext_sink = geompy.MakeTranslation(ext_sink , (n_active * base_x - ext_sink_x) / 2 ,-10 , -ext_sink_z)
    partition_array = np.append(partition_array , ext_sink)

    final_partition = geompy.MakePartition(partition_array.tolist(),  [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

    return(final_partition , ridge)


def create_cartridge(sc_data, geompy): #create the cartridge block
    [cart_x , cart_y , cart_z] = sc_data.cartridge_dim * 1000
    cartridge = geompy.MakeBoxDXDYDZ(cart_x, cart_y, cart_z)
    return(cartridge)

def create_chuck(sc_data , geompy): #create the chuck with walls to hold the cartridge

    [cart_x, cart_y, cart_z] = sc_data.cartridge_dim

    chuck_array = []
    chuck_base = geompy.MakeBoxDXDYDZ(300*10**3 , 100*10**3, 20*10**3)
    chuck_array = np.append(chuck_array , chuck_base)
    chuck_wall = geompy.MakeBoxDXDYDZ(300*10**3 , 2*10**3, 5*10**3)

    chuck_wall1 = geompy.MakeTranslation(chuck_wall, 0 , 0 , 20*10**3)
    chuck_wall2 = geompy.MakeTranslation(chuck_wall, 0 , (2 + cart_y)*10**3 , 20*10**3)
    chuck_array = np.append(chuck_array , chuck_wall1)
    chuck_array = np.append(chuck_array , chuck_wall2)
    chuck = geompy.MakePartition(chuck_array.tolist(), [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
    return(chuck)
    

def create_mesh(sc_data , fin_partition , sub_mesh_group , smesh):

    fine_mesh = sc_data.r_mesh ; coarse_mesh = sc_data.bdy_mesh

    Mesh_1 = smesh.Mesh(fin_partition,'Mesh_1')
    GMSH = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH)
    Gmsh_Parameters = GMSH.Parameters()
    Gmsh_Parameters.Set2DAlgo( 0 )
    Gmsh_Parameters.SetMinSize( 0 ) #minimum mesh size
    Gmsh_Parameters.SetMaxSize( 1e+22 ) #maximum mesh size 
    Gmsh_Parameters.SetSizeFactor( fine_mesh) #RIDGE MESH
    Gmsh_Parameters.SetIs2d( 0 )
    #Auto_group_for_Sub_mesh_1_2 = Mesh_1.GroupOnGeom(Auto_group_for_Sub_mesh_1,'Auto_group_for_Sub-mesh_1',SMESH.VOLUME)
    GMSH_1_1 = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH , geom = sub_mesh_group)
    Gmsh_Parameters_1 = GMSH_1_1.Parameters()
    Gmsh_Parameters_1.Set2DAlgo( 0 )
    Gmsh_Parameters_1.SetMinSize( 0 ) #minimum mesh size
    Gmsh_Parameters_1.SetMaxSize( 1e+22 ) #maximum mesh size
    Gmsh_Parameters_1.SetSizeFactor( coarse_mesh ) #BODY MESH
    Gmsh_Parameters_1.SetIs2d( 0 )
    
    isDone = Mesh_1.Compute()
    return(Mesh_1)

def find_face_of_god(smesh , group_array):

    for i  in range(0 , len(group_array)):
            BBox = smesh.GetBoundingBox(group_array[i])
            temp_z = 1
            if BBox.minZ == BBox.maxZ and BBox.minZ < temp_z:
                temp_z = BBox.minZ
                f_o_g = i + 1
    return(int(f_o_g))

def create_solid_array(Mesh_1 , partition_exploded):

    solid_array = []
    for i in range( 0 , len(partition_exploded)):
        solid_array = np.append(solid_array , Mesh_1.GroupOnGeom(partition_exploded[i],'Solid_0',SMESH.VOLUME)) 
    solid_array =  Mesh_1.GetGroups()
    group_array = Mesh_1.GetMesh().FaceGroupsSeparatedByEdges( 30, 0, 0 ) #seperated by edges for the mesh

    return(solid_array , group_array)

def new_mesh_ext_sink(data, arg0 ): # (ridge mesh , body mesh)

    geompy = geomBuilder.New()
    O = geompy.MakeVertex(0, 0, 0)
    OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
    OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
    OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

    #multi_ridge = create_ridges (data , geompy)
    chip, multi_ridge = create_device( data ,  geompy , back_facet_trench=False)

    if data.cartridge_mat!=0:
        adams_partition = []
        for i in range(0,data.n_chips):
            x_pos = (63 + i*173/2)*1000
            z_pos = (20+12)*1000 + data.ext_sink_dim[2]
            chip1 = geompy.MakeTranslation(chip, x_pos , 2.1*1000, z_pos )
            adams_partition = np.append(adams_partition , chip1)
        cartridge = create_cartridge(data, geompy)
        cartridge = geompy.MakeTranslation(cartridge, (300 - 180)*1000/2, 2*1000, 20*1000)
        adams_partition = np.append(adams_partition , cartridge)
        chuck = create_chuck(data, geompy)
        adams_partition = np.append(adams_partition , chuck)
        adams_partition = geompy.MakePartition(adams_partition.tolist(),  [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
    else:
        adams_partition = chip

    partition_exploded = geompy.ExtractShapes(adams_partition, geompy.ShapeType["SOLID"], False) #exploded the object into a large array
    sub_mesh_auto_group = geompy.CreateGroup(adams_partition, geompy.ShapeType["SOLID"])
    if data.cartridge_mat!=0:
        geompy.UnionList(sub_mesh_auto_group, partition_exploded[-4:]) #reverse order goes: Sink, base , ridges - chronological order
    else:
        geompy.UnionList(sub_mesh_auto_group, partition_exploded[-1:]) #reverse order goes: Sink, base , ridges - chronological order

    geompy.addToStudy( O, 'O' )
    geompy.addToStudy( OX, 'OX' )
    geompy.addToStudy( OY, 'OY' )
    geompy.addToStudy( OZ, 'OZ' )
    geompy.addToStudy( multi_ridge, 'multi ridge' )
    geompy.addToStudy(chip , 'final partition')

    for i in range( 0 , len(partition_exploded)):
        geompy.addToStudyInFather( chip, partition_exploded[i], f'Solid_{i}' )
    
    smesh = smeshBuilder.New()
    Mesh_1 = create_mesh(data , adams_partition , sub_mesh_auto_group, smesh)
    solid_array , group_array = create_solid_array(Mesh_1 , partition_exploded)

    f_o_g = find_face_of_god(smesh , group_array)   

    try:
      Mesh_1.ExportUNV( f'C:/ElmerFEM/ElmerFEM/bin/{arg0}/temp_save.unv', 0 )
      pass
    except:
      print('ExportUNV() failed. Invalid file name?')

    write_boundary_conds(f_o_g, arg0,  T_sink = device.T_sink) #can add convection boundaries if interested

    return(f_o_g)
    

arg1 = sys.argv[1] #project name - can also take in more arguments if necessary

V = float(sys.argv[2]) #input variable
sweeping_V = int(sys.argv[3])



pandas_data = pd.read_csv('C:/Projects/Perry_run/input_csv.csv').to_numpy()
device = MySemiconductor(pandas_data) #put data into MySemiconductor class
device.device_arb_parameter = -1 #placeholder to have the device ridge slightly not off the edge of the whole thing

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
elif sweeping_V == 9:
    device.device_arb_parameter = V

new_mesh_ext_sink(device, arg1)

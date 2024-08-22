#script for multiple 2mm cavity ridges like before
#from InputFolder.MySemiconductorClass import MySemiconductor
import sys
import salome
import numpy as np
from datetime import datetime
import pandas as pd
import os.path



sys.path.append('C:/Projects/bin')
from origin_write import MySemiconductor


#from Perry_active import device_hold


salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0,r'C:/Users/pmillar/Documents/Internship 2024/Week 4/Salome')


import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

def write_boundary_conds(boundary1 ,path, T_sink = 25, boundaries2 = 0): #heat sink temperature - though is not necessarily that complicated to change the input into an arry
    #here we need the boundary number of the lowest boundary (which was far too much effort to determing from the stupid Salome simulations)
    bound_cond = f'\n\nBoundary Condition 1\n  Target Boundaries(1) = {boundary1}\n  Name = "Heat Sink"\n  Temperature = {T_sink}\nEnd'
    if boundaries2 != 0:
        bound_cond2 = f'\n Boundary Condition 2 \n  Target Boundaries(1) = {boundaries2} \n Name = "Convection" \n External Temperature = 25 \n  Heat Transfer Coefficient = 25 \n End \n'
        bound_cond = bound_cond + boundaries2
    else:
        pass
    file = open(path , 'a') ; file.write(bound_cond) ; file.close()
    return()



def new_mesh_ext_sink(data , yhanger = 0): # (ridge mesh , body mesh)




    if True:
        #n_ridges ,y_hang = 200 , sinkz = 450 , sinkx = 3500 , sinky = 2400 , m1 = 0.03 , m2 = 0.8
    
    
        geompy = geomBuilder.New()
        
        O = geompy.MakeVertex(0, 0, 0)
        OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
        OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
        OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
    
    
        
        n_active = data.n_ridges #number of active regions / ridges
        thicknesses = data.r_heights #thickness for each materials
        widths = data.r_widths
        lengths = data.r_lengths
        
        
        [box_x , box_y , box_z] = data.device_dim #device dimensions in um
    
        ridge_mesh_fac = data.r_mesh
        body_mesh_fac = data.bdy_mesh
    
        [sinkx, sinky , sinkz] = data.ext_sink_dim
    
    
        
        
        if thicknesses.size != widths.size: #again shitty little if statement to make sure every ridge has as a corresponding width
            blow_up = 1/0
        else:
            pass
        
        z_active = 0 #z_height position of active position, should probably have a check that this region is in fact within the device size
        
        
        heater_y = 2000 #um length of the heater on top of the device
        #again assume the heater is the upper most layer on the ridge and that it is centered on the ridge. 
        
        
        z_ridge = data.z_ridge # z-position of the bottom point of the ridge
        
        #get the device size check done now to make the loop below a bit cleaner
        #need to check that x-span of all the active regions is not greater than the box dimensions
        #dividing by zero just stops the code with an error, can't be bothered trying to import stuff to properly exit - which is essentially the same thing
        
        
        
        
        
        
        partition_dummy = []
        
        
        #create ridges and add them to an array, inlcuding the vertical shift to layer them on each other
        #same process unless different ridges have different layers / thicknesses
        
        
        z_height = 0 #initialize height for ridges to sit on
        for i in range(0,len(thicknesses)):
            if i ==0: #first ridge
                ridge = geompy.MakeBoxDXDYDZ(widths[i], lengths[i] ,  thicknesses[i]) #get ridge depth AND height
                ridge = geompy.MakeTranslation(ridge, -widths[i]/2,(box_y - lengths[i]) / 2, 0) #centre the origin 
                partition_dummy = np.append(partition_dummy ,  ridge)#
                geompy.addToStudy(ridge , f'ridge{i}')
            elif i == len(thicknesses)-1: #last ridge ie heater ridge - centre on the device
                ridge = geompy.MakeBoxDXDYDZ(widths[i], lengths[i] , thicknesses[i]) #create box
                ridge = geompy.MakeTranslation(ridge, -widths[i]/2, (box_y - lengths[i]) / 2, z_height + thicknesses[i-1]) #shift all half width back to center origin at symmetry centre
                partition_dummy = np.append(partition_dummy ,  ridge)
                z_height += thicknesses[i-1]
                geompy.addToStudy(ridge , f'ridge{i}')  
            else:
                ridge = geompy.MakeBoxDXDYDZ(widths[i], box_y , thicknesses[i]) #create box
                ridge = geompy.MakeTranslation(ridge, -widths[i]/2,0, z_height + thicknesses[i-1]) #shift all half width back to center origin at symmetry centre
                partition_dummy = np.append(partition_dummy ,  ridge)
                z_height += thicknesses[i-1]
                geompy.addToStudy(ridge , f'ridge{i}')
        
        partition_dummy = partition_dummy.tolist()
        
        multi_ridge = geompy.MakePartition(partition_dummy, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0) #ridge - final product for us now
        
        
        Base = geompy.MakeBoxDXDYDZ(box_x * n_active , box_y, box_z) # create the the base

        #300um Falcon base
        
    
        
        #now need to cut the ridge out of the box
        
        partition_dummy2 = []
        
        for i in range(0, n_active):
            x_pos = 125 + (50 * i) #want each active region in its own device ~250um so separation is 250um in total
            multi_ridge_cut = geompy.MakeTranslation(multi_ridge , x_pos, 0 , z_ridge) #this is actually our final object we want to work with
            Base = geompy.MakeCut(Base, multi_ridge_cut, True) #Final box is sorted now
            partition_dummy2 = np.append(partition_dummy2 , multi_ridge_cut) #add ridge in new location to the partition
        
        
        partition_dummy2 = np.append(partition_dummy2 , Base)

        no_in_submesh = 1 #need a counter for how many objects are going to be in the device

        back_facet_array = np.array([0.219 ,0.097 , 0.165,0.097 , 0.165,0.097 , 0.165,0.097]) #back facet thickness
        
        if data.ffront_mat != 0 or data.fback_mat !=0:#no facet coating present
            #create and move front facets
            front_facet = geompy.MakeBoxDXDYDZ(box_x * n_active , data.ffront_y, box_z) #create front facet to cover front of the box
            front_facet = geompy.MakeTranslation(front_facet ,0 ,-data.ffront_y , 0 )#move infront of the device
            partition_dummy2 = np.append(partition_dummy2 , front_facet )
            #create and move back facets


            for bf in range(0,len(back_facet_array)):
                if bf==0:
                    back_facet = geompy.MakeBoxDXDYDZ(box_x * n_active , back_facet_array[bf], box_z) #create back facet to cover front of the box
                    back_facet = geompy.MakeTranslation(back_facet ,0 ,box_y , 0 )#move infront of the device
                    partition_dummy2 = np.append(partition_dummy2 , back_facet)
                else:
                    back_facet = geompy.MakeBoxDXDYDZ(box_x * n_active , back_facet_array[bf], box_z) #create back facet to cover front of the box
                    back_facet = geompy.MakeTranslation(back_facet , 0 , box_y + np.sum(back_facet_array[:bf]) , 0 )#move infront of the device
                    partition_dummy2 = np.append(partition_dummy2 , back_facet)


            no_in_submesh += len(back_facet_array) + 1 #tells the submesh that we also have two other objects to include
        else:
            pass


        if device.ext_sink_mat != 0: #create external heat sink
            AlN_base = geompy.MakeBoxDXDYDZ(sinkx, sinky, sinkz)
            AlN_base = geompy.MakeTranslation(AlN_base , -0.5* (sinkx - box_x*n_active),-0.5*(sinky - box_y),-sinkz) #subtract the hanging distance off
            partition_dummy2 = np.append(partition_dummy2 , AlN_base)
            no_in_submesh += 1 #adds another object (if so) to be sub-meshed
            #-0.5*(sinky - box_y)
        else:
            pass
    
        
        
        
        final_partition = geompy.MakePartition(partition_dummy2.tolist(),  [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
        
        partition_exploded = geompy.ExtractShapes(final_partition, geompy.ShapeType["SOLID"], False) #exploded the object into a large array
        #partition_exploded = geompy.GetExistingSubObjects(final_partition, False)

        sub_mesh_auto_group = geompy.CreateGroup(final_partition, geompy.ShapeType["SOLID"])
        

        #everything from partition exploded that goes into the line below is included in the submesh
        #we want anything that is NOT a ridge to be here ie facets, base , heat sinks 
        #sub mesh is actually coarser not finer, but we need the submesh to be calculated first 
        #to avoid overmeshing on the base and at the facets
        geompy.UnionList(sub_mesh_auto_group, partition_exploded[-no_in_submesh:])
    
        
        #partition together just the ridge. Assuming the ridge has the same layer geometry for each ridge
        #[a,b,c,d] =  geompy.ExtractShapes(final_partition, geompy.ShapeType["SOLID"], True)
        #[a,b,c,d] =  geompy.GetExistingSubObjects(final_partition, False)
        
        geompy.addToStudy( O, 'O' )
        geompy.addToStudy( OX, 'OX' )
        geompy.addToStudy( OY, 'OY' )
        geompy.addToStudy( OZ, 'OZ' )
        
        geompy.addToStudy( multi_ridge, 'multi ridge' )
        geompy.addToStudy(final_partition , 'final partition')
        
        for i in range( 0 , len(partition_exploded)):
            geompy.addToStudyInFather( final_partition, partition_exploded[i], f'Solid_{i}' )
        
        
        
        
        #geompy.addToStudyInFather(final_partition , a, 'a')
        #geompy.addToStudyInFather(final_partition , b, 'b')
        #geompy.addToStudyInFather(final_partition , c, 'c')
        #geompy.addToStudyInFather(final_partition , d, 'd')
    
    
        
        
        #end of script before meshing
        #a,b,c,d,e....should correspond to the total number of objects. ie if there were ridges, r each with layers, l the number of objects (letters) would need to be r*l+1 (+1 for the base)
        import  SMESH, SALOMEDS
        from salome.smesh import smeshBuilder
        
        
        smesh = smeshBuilder.New()
#sme    sh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                     # multiples meshes built in parallel, complex and numerous mesh edition (performance)
    
    
    
        Mesh_1 = smesh.Mesh(final_partition,'Mesh_1')
        GMSH = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH)
        Gmsh_Parameters = GMSH.Parameters()
        Gmsh_Parameters.Set2DAlgo( 0 )
        Gmsh_Parameters.SetMinSize( 0 ) #minimum mesh size
        Gmsh_Parameters.SetMaxSize( 1e+22 ) #maximum mesh size 
        Gmsh_Parameters.SetSizeFactor( ridge_mesh_fac) #RIDGE MESH
        Gmsh_Parameters.SetIs2d( 0 )
    
    
        #Auto_group_for_Sub_mesh_1_2 = Mesh_1.GroupOnGeom(Auto_group_for_Sub_mesh_1,'Auto_group_for_Sub-mesh_1',SMESH.VOLUME)
        
        GMSH_1_1 = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH , geom = sub_mesh_auto_group)
        Gmsh_Parameters_1 = GMSH_1_1.Parameters()
        Gmsh_Parameters_1.Set2DAlgo( 0 )
        Gmsh_Parameters_1.SetMinSize( 0 ) #minimum mesh size
        Gmsh_Parameters_1.SetMaxSize( 1e+22 ) #maximum mesh size
        Gmsh_Parameters_1.SetSizeFactor( body_mesh_fac ) #BODY MESH
        Gmsh_Parameters_1.SetIs2d( 0 )
        
        isDone = Mesh_1.Compute()

        #XXXSubmesh = GMSH_1_1.GetSubMesh()
        
        solid_array = []
        
        for i in range( 0 , len(partition_exploded)):
            solid_array = np.append(solid_array , Mesh_1.GroupOnGeom(partition_exploded[i],'Solid_0',SMESH.VOLUME)) 
        
        Sub_mesh_1 = GMSH_1_1.GetSubMesh()
        
        solid_array =  Mesh_1.GetGroups()
    
        ##get shape bounding box info - criterion for the lower most face
        
        group_array = Mesh_1.GetMesh().FaceGroupsSeparatedByEdges( 30, 0, 0 ) #seperated by edges for the mesh
    
    
        face_ind = np.array([])
        face_z = np.array([])

        open_boundary = np.array([])
    
        for i  in range(0 , len(group_array)):
            BBox = smesh.GetBoundingBox(group_array[i])

            if BBox.minZ == BBox.maxZ and (BBox.maxY - BBox.minY) > box_y/2:
                #face_ind = np.append(face_ind, BBox.minZ, axis = 0)
                face_ind = np.append(face_ind, i-1) ; face_z = np.append(face_z , BBox.minZ)
            else:
                pass


            if data.ffront_mat != 0 or data.fback_mat !=0: #facets are present - go on to include convection boundary conditions

                if BBox.maxZ == BBox.minZ and BBox.maxZ == box_z and BBox.maxZ == box_y: #top faces of the device
                    open_boundary = np.append(open_boundary , i)
                else:
                    pass
            
                if BBox.maxX == BBox.minX and BBox.minX == 0: #get x+ and x- open faces
                    pass
                    #open_boundary = np.append(open_boundary , i)
                elif BBox.maxX == BBox.minX and BBox.maxX == n_active*box_x:
                    #open_boundary = np.append(open_boundary , i)
                    pass
                if BBox.maxY == BBox.minY and BBox.maxY == box_y + data.fback_y: #open back facet
                    #open_boundary = np.append(open_boundary , i)
                    pass
                elif  BBox.maxY == BBox.minY and BBox.minY == -data.ffront_y: #open front facet
                    #open_boundary = np.append(open_boundary , i)
                    pass
                else:
                    pass




        face_of_god = int(face_ind[np.argmin(face_z)] + 2)
    else:
        pass
    

    if data.ffront_mat != 0 or data.fback_mat !=0:
        open_string = ' '
        for x in open_boundary:
            open_string += str(int(x + 1)) + ' '
    else:
        open_string = 0

    try:
      Mesh_1.ExportUNV( f'C:/ElmerFEM/ElmerFEM/bin/{arg0}/temp_save.unv', 0 )
      #Mesh_1.ExportUNV( f'C:/ElmerFEM/ElmerFEM/bin/{arg0}/{open_string}.unv', 0 )
      pass
    except:
      print('ExportUNV() failed. Invalid file name?')

      

    mypath = f'C:/ElmerFEM/ElmerFEM/bin/{arg0}/'
    sif_directory = os.path.join(mypath , 'case.sif')
    write_boundary_conds(face_of_god, sif_directory , T_sink = device.T_sink, boundaries2= 0) #can add convection boundaries if interested

    return(face_of_god)
    

#C:\Users\pmillar\Documents\Internship 2024\Week 8
#C:\SALOME-9.12.0\W64\Python\InputFolder

arg0 = sys.argv[1] #project name - can also take in more arguments if necessary

V = float(sys.argv[2]) #input variable
sweeping_V = int(sys.argv[3])

pandas_data = pd.read_csv('C:/Projects/bin/input_csv.csv').to_numpy()
device = MySemiconductor(pandas_data) #put data into MySemiconductor class

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

new_mesh_ext_sink(device)

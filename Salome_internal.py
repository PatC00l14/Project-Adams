import sys
import salome
import numpy as np
import pandas as pd
import os.path

sys.path.append('C:/Projects/Perry_run')
from origin_write import *


salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

import GEOM
from salome.geom import geomBuilder 
import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

def write_ridge_bodies( device, project_name, scc = 0): 
    body_string = '\n'

    n_r = device.n_ridges #number of ridges per chip
    n_l = device.n_layers #number of layers per ridge

    if device.current_model == 0:
        for r in range(0, n_r): #ridge index
            count = 1
            for l in range(1,n_l+1): #layer index
                num = r*(n_l) + l
                if device.r_heat_power[l-1] *device.r_onoff[r] !=0: #this means there IS heating in this layer
                    body_string = body_string + write_ind_body(num, device.r_materials[l-1] ,count) + f'\n'
                    count += 1
                else:
                    body_string = body_string + write_ind_body(num, device.r_materials[l-1] , 0) + f'\n'
        num = (n_r*n_l) + 1

    elif device.current_model == 1:
        for r in range(0, n_r): #ridge index
                count = 1
                for l in range(1,n_l+1): #layer index
                    num = r*(n_l) + l
                    body_string = body_string + write_ind_body(num, device.r_materials[l-1] , 0) + f'\n'   
    else:
        raise('Error: wrong input type of SCC on (1) or off (0)')

    file = open(f'C:/ElmerFEM/ElmerFEM/bin/{project_name}/case.sif' , 'a')
    file.write(f'{body_string}')
    file.close()
    return()

def write_ind_bodyx(ind , material , body_force, project_name, name = 'Body Property'): # write an individual body force term
 
    material = int(material)

    if body_force ==0:
        bdy_string00 = f'Body {ind} \n  Target Bodies(1) = {ind}  \n  Name = "{name} {ind}" \n  Equation = 1 \n  Material = {material} \n  Initial condition = 1 \nEnd \n'
    else:
        bdy_string00 = f'Body {ind} \n  Target Bodies(1) = {ind}  \n  Name = "{name} {ind}" \n  Equation = 1 \n  Material = {material} \n  Body Force = {body_force} \n  Initial condition = 1 \nEnd \n'
    my_file = open(f'C:/ElmerFEM/ElmerFEM/bin/{project_name}/case.sif' , 'a')
    my_file.write(bdy_string00)
    my_file.close()
    return (bdy_string00)



def write_boundary_conds(boundary1 , arg0 ,dat, bound_scc = 0, potential = 1.7, ro_i = (0.75*10**7)): #heat sink temperature - though is not necessarily that complicated to change the input into an arry
    #here we need the boundary number of the lowest boundary (which was far too much effort to determing from the stupid Salome simulations)
    mypath = f'C:/ElmerFEM/ElmerFEM/bin/{arg0}/'
    path = os.path.join(mypath, 'case.sif')
    
    if dat.current_model ==0:
        bound_cond1 = f'\n\nBoundary Condition 1\n  Target Boundaries(1) = {boundary1}\n  Name = "Heat Sink"\n  Temperature = {dat.T_sink}\nEnd'
    
    elif dat.current_model ==1:
        bound_nums = int(len(bound_scc))
        scc_string = np.array2string(bound_scc)[1:-1]
        bound_cond1 = f'\n\nBoundary Condition 1\n  Target Boundaries({bound_nums}) = {scc_string}\n  Name = "Potential Bondary"\n  Potential = {potential}\n  Current Density =  {ro_i}\nEnd'        
        #bound_cond1 = f'\n\nBoundary Condition 1\n  Target Boundaries({bound_nums}) = {scc_string}\n  Name = "Potential Bondary"\n  Potential = {potential}\nEnd'
        bound_cond1 = bound_cond1 + f'\n\nBoundary Condition 2\n  Target Boundaries(1) = {boundary1}\n  Name = "GROUND"\n  Potential = 0\nEnd'
    else:
        raise Exception('Error: wrong input type of SCC on (1) or off (0)')
    
    try:
        file = open(path , 'a') 
        file.write(bound_cond1)
        file.close()
    except:
        raise Exception("Error: unable to write boundary conditions - perhaps face index was not found")
 
    
    return()


def create_ridges(sc_data , geompy, middle_y = True): #this function creates the combined ridge structure
    
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

        if middle_y:
            ridge = geompy.MakeTranslation(ridge , -widths[i]/2 , (base_y - lengths[i])/2 , np.sum(thicknesses[:i]))
        else:
            ridge = geompy.MakeTranslation(ridge , -widths[i]/2 , 0 , np.sum(thicknesses[:i]))
        geompy.addToStudy(ridge , f'ridge{i}')
        ridge_partition = np.append(ridge_partition , ridge)

    combined_ridge = geompy.MakePartition(ridge_partition.tolist(), [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

    #need to find the ethch stop location
    
    return(combined_ridge)


def write_partition_XXDEADDONTUSE(obj, part_arr, count,  mat, proj_name, name = 'Body Property', bdy = 0, repeat = 1):
    #function to add to the partition and write to the case.sif file at the same time to keep track of the order of materials added
    for i in range(repeat):
        part_arr = np.append(part_arr, obj)
        write_ind_bodyx(count, mat, bdy, proj_name, name = name)
        count += 1
    return(part_arr, count)

def write_partition(count,  mat, proj_name, name = 'Body Property', bdy = 0, repeat = 1):
    #function to add to the partition and write to the case.sif file at the same time to keep track of the order of materials added
    for i in range(repeat):
        write_ind_bodyx(count, mat, bdy, proj_name, name = name)
        count += 1
    return(count)

def get_wirebond_index(sc_data):
    #I-A-I-I-A
    # au - [6,9] ; ins - [5,7,8] --n = 1
    # au - [n*5+1,n*5+4] ; ins - [n*5,n*5+2,N*5+3]  -- n = n
    n_active = sc_data.n_ridges
    n0 = sc_data.n_layers * n_active
    
    if sc_data.bfm !=0:
        au_inds = [n0+1, n0+4]; insul_inds = [n0,n0+2,n0+3] 
        if  n_active > 1:
            for i in range(1,n_active):
                n = n0 + i*5
                extendAU = [n+0, n+3] ; extendINS = [n+1, n+2, n+4]
                au_inds.extend(extendAU) ; insul_inds.extend(extendINS)

    elif sc_data.bfm ==0:
        au_inds = [n0+1] ; insul_inds = [n0,n0+2]
        if  n_active > 1:
            for i in range(1,n_active):
                n = n0 + i*3
                extendAU = [n] ; extendINS = [n+1, n+2]
                au_inds.extend(extendAU) ; insul_inds.extend(extendINS)
    return(au_inds, insul_inds)

def get_current_model_wirebond_index(sc_data):

    n_active = sc_data.n_ridges
    n0 = sc_data.n_layers * n_active
    if sc_data.bfm !=0:
        au_inds = (np.arange(n0, n0+2*n_active,1).astype(int)).tolist()
    elif sc_data.bfm ==0:
        au_inds = (np.arange(n0, n0+1*n_active,1).astype(int)).tolist()

    return(au_inds)

def cut_bfm_trenches(data, geompy, obj_to_cut):

    [base_x , base_y , base_z] = data.device_dim
    [trench_x, trench_y , trench_z] = data.trench_dim

    bfm_trench1 = geompy.MakeBoxDXDYDZ(base_x*(data.n_ridges+1), trench_x, trench_z)
    bfm_trench1 = geompy.MakeTranslation(bfm_trench1, -base_x/2, base_y, base_z - trench_z)
    bfm_trench2 =  geompy.MakeTranslation(bfm_trench1, -base_x/2, base_y -20 , base_z - trench_z)    
    obj_to_cut = geompy.MakeCut(obj_to_cut, bfm_trench1)
    obj_to_cut = geompy.MakeCut(obj_to_cut, bfm_trench2)
    
    return(obj_to_cut)

def create_insulating_layer(data, geompy, base, base0, e_stop, pro_name , count = 0, write = True):

    empty_part = np.array([])
     
    [base_x , base_y , base_z] = data.device_dim
    [trench_x, trench_y , trench_z] = data.trench_dim
    n_ridges = data.n_ridges

    h_in = 0.3 #insulation thickness - will be assigned at some point
    h_au = 0.3 #Gold thickness - can be assigned now but I cba at the minute
    r_tx = 22 #ridge trench width
    r_wx = np.min(data.r_widths) #width of the ridge
    r_space = 0.5 * (r_tx - r_wx) #width of space from ridge edge to trench edge
    #                           _______
    #       --------|          |-r_wx-|          |--------
    #               |_r_space_|       |_________|
    #                           r_tx            

    insul_layer = geompy.MakeBoxDXDYDZ(base_x * (data.n_ridges), base_y+data.bfm, base_z + h_in)
    insul_layer = geompy.MakeCut(insul_layer, base0)    
    insul_layer = geompy.MakeCut(insul_layer, base)
    
    insul_trench = geompy.MakeBoxDXDYDZ(r_space - 2*h_in, base_y, e_stop)
    insul_trench = geompy.MakeTranslation(insul_trench, - 0.5 * (r_space - 2*h_in), 0, base_z - e_stop + h_in)

    au_layer = geompy.MakeBoxDXDYDZ(base_x* (data.n_ridges), base_y+data.bfm, base_z + h_in + h_au)
    au_layer = geompy.MakeCut(au_layer, base0)
    au_layer = geompy.MakeCut(au_layer, base)

    au_trench = geompy.MakeBoxDXDYDZ(r_space - 2*(h_in+h_au), base_y, e_stop)
    au_trench = geompy.MakeTranslation(au_trench, - 0.5 * (r_space - 2*(h_in+h_au)), 0, base_z - e_stop + h_in + h_au)
    
    for i in range(data.n_ridges):
        #cut right side of ridge
        x_pos1 = base_x * ( i + 0.5) + 0.5 * (r_wx + r_space)
        ins_cut = geompy.MakeTranslation(insul_trench, x_pos1, 0, 0)
        insul_layer = geompy.MakeCut(insul_layer, ins_cut)
        au_cut = geompy.MakeTranslation(au_trench,  x_pos1, 0 , 0 )
        au_layer = geompy.MakeCut(au_layer, au_cut)

        #cut left side of ridge
        x_pos2 = base_x * ( i + 0.5) - 0.5 * (r_wx + r_space)
        ins_cut = geompy.MakeTranslation(insul_trench, x_pos2, 0, 0)
        insul_layer = geompy.MakeCut(insul_layer, ins_cut)
        au_cut = geompy.MakeTranslation(au_trench,  x_pos2, 0 , 0 )
        au_layer = geompy.MakeCut(au_layer , au_cut)

        #cut top of ridge for p-type metal to sit 
        x_pos3 = base_x * ( i + 0.5)
        ins_cut = geompy.MakeTranslation(insul_trench, x_pos3, 0, -h_in+e_stop)
        insul_layer = geompy.MakeCut(insul_layer, ins_cut)
    
    if write:
        if data.current_model == 0:
            au_inds, insul_inds = get_wirebond_index(data)
            for au in au_inds:
                count = write_partition(au, 2, pro_name, name = 'Au cap')
            for ins in insul_inds:
                count = write_partition(ins, 13, pro_name, name = 'insul cap')
    
        elif data.current_model ==1:
            au_inds = get_current_model_wirebond_index(data)
            for au in au_inds:
                count = write_partition(au, 2, pro_name, name = 'au cap')

    return(insul_layer, au_layer, count)


def create_insulating_layer2_devopsstuff(data, geompy, base, base0, e_stop, pro_name , partition, count = 0, write = True):

    empty_part = np.array([])
     
    [base_x , base_y , base_z] = data.device_dim
    [trench_x, trench_y , trench_z] = data.trench_dim
    n_ridges = data.n_ridges

    h_in = data.insul_z #insulation thickness - will be assigned at some point
    h_au = data.au_cap #Gold thickness - can be assigned now but I cba at the minute
    r_tx = 22 #ridge trench width
    r_wx = np.min(data.r_widths) #width of the ridge
    r_space = 0.5 * (r_tx - r_wx) #width of space from ridge edge to trench edge
    #                           _______
    #       --------|          |-r_wx-|          |--------
    #               |_r_space_|       |_________|
    #                           r_tx            

    insul_layer = geompy.MakeBoxDXDYDZ(base_x -(trench_x + 20), base_y, base_z + h_in)
    insul_layer = geompy.MakeTranslation(insul_layer, 0.5 * (trench_x + 20), 0, 0)
    insul_layer = geompy.MakeCut(insul_layer, base0)    
    insul_layer = geompy.MakeCut(insul_layer, base)
    insul_trench = geompy.MakeBoxDXDYDZ(r_space - 2*h_in, base_y, e_stop)
    insul_trench = geompy.MakeTranslation(insul_trench, - 0.5 * (r_space - 2*h_in), 0, base_z - e_stop + h_in)

    au_layer = geompy.MakeBoxDXDYDZ(base_x-(trench_x + 20), base_y, base_z + h_in + h_au)
    au_layer = geompy.MakeTranslation(au_layer,  0.5 * (trench_x + 20), 0, 0)
    au_layer = geompy.MakeCut(au_layer, base0)
    au_layer = geompy.MakeCut(au_layer, base)
    au_trench = geompy.MakeBoxDXDYDZ(r_space - 2*(h_in+h_au), base_y, e_stop)
    au_trench = geompy.MakeTranslation(au_trench, - 0.5 * (r_space - 2*(h_in+h_au)), 0, base_z - e_stop + h_in + h_au)
    
    i = 0
    #cut right side of ridge
    x_pos1 = base_x * ( i + 0.5) + 0.5 * (r_wx + r_space)
    ins_cut = geompy.MakeTranslation(insul_trench, x_pos1, 0, 0)
    insul_layer = geompy.MakeCut(insul_layer, ins_cut)
    au_cut = geompy.MakeTranslation(au_trench,  x_pos1, 0 , 0 )
    au_layer = geompy.MakeCut(au_layer, au_cut)
    #cut left side of ridge
    x_pos2 = base_x * ( i + 0.5) - 0.5 * (r_wx + r_space)
    ins_cut = geompy.MakeTranslation(insul_trench, x_pos2, 0, 0)
    insul_layer = geompy.MakeCut(insul_layer, ins_cut)
    au_cut = geompy.MakeTranslation(au_trench,  x_pos2, 0 , 0 )
    au_layer = geompy.MakeCut(au_layer , au_cut)
    #cut top of ridge for p-type metal to sit 
    x_pos3 = base_x * ( i + 0.5)
    ins_cut = geompy.MakeTranslation(insul_trench, x_pos3, 0, -h_in+e_stop)
    insul_layer = geompy.MakeCut(insul_layer, ins_cut)

    if data.bfm !=0:
        if 1==2:
            au_layer = cut_bfm_trenches(data, geompy, au_layer)
            insul_layer = cut_bfm_trenches(data, geompy, insul_layer)
        else:
            bfm_trench1 = geompy.MakeBoxDXDYDZ(base_x*(data.n_ridges+1), trench_x, trench_z)
            bfm_trench1 = geompy.MakeTranslation(bfm_trench1, -base_x/2, base_y, base_z - trench_z)
            bfm_trench2 =  geompy.MakeTranslation(bfm_trench1, -base_x/2, base_y -20 , base_z - trench_z)

            insul_layer = geompy.MakeCut(insul_layer, bfm_trench1)
            au_layer = geompy.MakeCut(au_layer, bfm_trench1)

            insul_layer = geompy.MakeCut(insul_layer, bfm_trench2)
            au_layer = geompy.MakeCut(au_layer, bfm_trench2)

    empty = np.array([])
    
    if write:
        for i in range(n_ridges):
            au_layer_update = geompy.MakeTranslation(au_layer, base_x*i, 0, 0)
            partition = geompy.append(partition, au_layer_update)
            count = write_partition(count, 2, pro_name, name = 'Au layer')
            insul_layer_update = geompy.MakeTranslation(insul_layer, base_x*i, 0, 0)
            partition = np.append(partition, insul_layer_update)
            count = write_partition(count, 13, pro_name, name = 'Insul layer')
        return(partition, count, insul_layer, au_layer )
    else:
        return(insul_layer, au_layer)


def create_wirebonds(data, geompy, partition, pro_name, n_wirebonds = 1, count = 0):
    
    wire_bond_partition = np.array([])
    [base_x, base_y, base_z] = data.device_dim
    if data.au_cap_mat ==0 or data.insul_mat ==0:
        return(0,0)
    else:
        w_unit = 20
        wire_bond = geompy.MakeBoxDXDYDZ(w_unit,w_unit,w_unit)
        for i in range(data.n_ridges):
            for j in range(1,n_wirebonds+1):
                x_pos = (0.25 + i)*base_x - w_unit/2
                y_pos = j*base_y / (n_wirebonds+1)
                wire_bond_set = geompy.MakeTranslation(wire_bond, x_pos, y_pos, base_z + data.au_cap+data.insul_z)
                partition = np.append(partition, wire_bond_set)
        
        count = write_partition(count, data.au_cap_mat, pro_name, repeat=data.n_ridges * n_wirebonds, name='Wirebonds')          

        return(partition, count)



def create_pdown_layer(data, geompy, base, base0, estop,  pro_name, count, partition, write = True):
    #pside down creation
    [base_x, base_y, base_z] = data.device_dim
    [trench_x, trench_y, trench_z] = data.trench_dim

    insul_layer , au_ignore = create_insulating_layer2_devopsstuff(data, geompy, base, base0, estop, pro_name, partition, write = False)
    
    empty = np.array([])

    au_body = geompy.MakeBoxDXDYDZ(base_x - (trench_x + 10), base_y, estop + data.au_cap+ data.insul_z)
    au_body = geompy.MakeTranslation(au_body, (trench_x + 10)/2, 0, base_z - estop)

    au_body = geompy.MakeCut(au_body, insul_layer, True)
    au_body = geompy.MakeCut(au_body, base, True)

    for i in range(data.n_ridges):
        au_bodyCUT = geompy.MakeTranslation(au_body, i*base_x, 0, 0) 
        count = write_partition(count, 2, pro_name, name = 'Thick Au filling')
        partition = np.append(partition, au_bodyCUT)

        insul_update = geompy.MakeTranslation(insul_layer, i*base_x, 0, 0 )
        count = write_partition(count, 13, pro_name, name = 'insulation layer', repeat= 2)
        partition = np.append(partition, insul_update)

    opening_dim = [70, 730, 0.2]

    au_inbetweenlayer = geompy.MakeBoxDXDYDZ(base_x*data.n_ridges, base_y, opening_dim[2])
    au_inbetweenlayer = geompy.MakeTranslation(au_inbetweenlayer, 0, 0, base_z+ data.au_cap+ data.insul_z)

    opening = geompy.MakeBoxDXDYDZ(opening_dim[0]/2, opening_dim[1], opening_dim[2])
    opening = geompy.MakeTranslation(opening, -opening_dim[0]/4, (base_y - opening_dim[1])/2, base_z+ data.au_cap+ data.insul_z)

    for i in range(0, 2* data.n_ridges):
        xpos = 0.5 *(i+0.5)*base_x
        opening_cut = geompy.MakeTranslation(opening,xpos, 0, 0 )
        au_inbetweenlayer = geompy.MakeCut(au_inbetweenlayer, opening_cut, True)

    count = write_partition(count, 13, pro_name, name = 'Isolation')
    partition = np.append(partition, au_inbetweenlayer)

    return(partition, count)
    
def create_thick_au_pads(data, geompy, partition, pro_name, count = 0):
    
    wire_bond_partition = np.array([])
    [base_x, base_y, base_z] = data.device_dim
    if data.au_cap_mat ==0 or data.insul_mat ==0:
        return(0)
    else:
        au_pad_z = 1.0
        au_pad_x = 60.0
        wire_bond = geompy.MakeBoxDXDYDZ(au_pad_x,base_y,au_pad_z)
    
        for i in range(data.n_ridges):
            x_pos = (0.1 + i)*base_x 
            wire_bond_set = geompy.MakeTranslation(wire_bond, x_pos,0, base_z + data.au_cap+data.insul_z)
            partition = np.append(partition, wire_bond_set)
        
        bdy_num = np.sum(1 for x in data.r_heat_power if x!=0)
        count = write_partition(count, data.au_cap_mat, pro_name, repeat=data.n_ridges, bdy = bdy_num + 1, name= 'Thick Au Pad')          
    
        return(partition, count)
    


def create_heater(data, geompy, partition, pro_name, count=0):

    [base_x, base_y, base_z] = data.device_dim
    [heat_x, heat_y, heat_z] = [np.min(data.r_widths),300, 0.5]
    [ins_x, ins_y, ins_z] = [np.min(data.r_widths), base_y, 0.5]

    heater_height = base_z + data.au_cap + data.insul_z

    insul = geompy.MakeBoxDXDYDZ(ins_x, ins_y, ins_z)
    insul = geompy.MakeTranslation(insul, -ins_x/2, 0 , heater_height)

    heater = geompy.MakeBoxDXDYDZ(heat_x, heat_y, heat_z)   
    heater = geompy.MakeTranslation(heater, -heat_x/2, 0.5*(base_y - heat_y), heater_height+ins_z) 

    empty_part = []

    for i in range(data.n_ridges):
        xpos = (0.5+i)*base_x
        heater_place = geompy.MakeTranslation(heater, xpos, 0, 0)
        ins_place = geompy.MakeTranslation(insul, xpos, 0, 0)

        bdy_num = int(sum( 1 for x in data.r_heat_power if x!=0 )) + 1

        partition = np.append(partition, ins_place)
        count = write_partition(count, data.insul_mat, pro_name, name = 'Heater insulation')

        partition = np.append(partition, heater_place)
        count = write_partition(count, data.insul_mat,pro_name, bdy = bdy_num, name = 'Heater')

    return(partition, count)

def create_submount(data, geompy, partition, count, pro_name):

    [ext_sink_x,ext_sink_y,ext_sink_z] = data.ext_sink_dim
    n_active = data.n_ridges
    [base_x, base_y, base_z] = data.device_dim

    empty = np.array([])

    if data.pside_down == 0:
         ext_sink = geompy.MakeBoxDXDYDZ(ext_sink_x,ext_sink_y,ext_sink_z) 
         ext_sink = geompy.MakeTranslation(ext_sink , (n_active * base_x - ext_sink_x) / 2 ,-100 , -ext_sink_z)
         count = write_partition(count, 6,pro_name,  name = 'Submount')
         partition = np.append(partition, ext_sink )

    elif data.pside_down == 1:
        ext_sink1 = geompy.MakeBoxDXDYDZ(ext_sink_x,ext_sink_y,ext_sink_z - 40)
        ext_sink1 = geompy.MakeTranslation(ext_sink1 , (n_active * base_x - ext_sink_x) / 2 ,-100 , base_z+data.au_cap+data.insul_z+0.2+40)

        ext_sink2 = geompy.MakeBoxDXDYDZ(ext_sink_x,ext_sink_y,40)
        ext_sink2 = geompy.MakeTranslation(ext_sink2 , (n_active * base_x - ext_sink_x) / 2 ,-100 , base_z+data.au_cap+data.insul_z+0.2)

        
        count = write_partition(count, 6,pro_name,  name = 'Submount main')
        count = write_partition(count, 6, pro_name, name = 'Submount buffer')
        partition = np.append(partition, ext_sink2) ; partition = np.append(partition, ext_sink1)

    return(partition, count)


def find_estop(data):
    ##quick method to find the estop height to determine the ridge trench depth
    r_widths = data.r_widths
    for i in range(1 , len(r_widths)):
        if r_widths[i] < r_widths[i-1]:
            num = i
            break
    return(np.sum(data.r_heights[num:]))

def create_device ( data, geompy, back_facet_trench = False, middle_y = True, extra_width = True, pro_name = 'y'):
    #this will create the a single chip with a submount. The ridge is created within this function

    ridge = create_ridges(data , geompy, middle_y = middle_y)

    partition_array = []

    [base_x , base_y , base_z] = data.device_dim
    [trench_x, trench_y , trench_z] = data.trench_dim
    [ext_sink_x,ext_sink_y,ext_sink_z] = data.ext_sink_dim
    n_active = data.n_ridges

    z_pos = base_z  - np.sum(data.r_heights)
    e_stop = find_estop(data)

    trench = geompy.MakeBoxDXDYDZ(trench_x , trench_y + data.bfm , trench_z+20)
    trench = geompy.MakeTranslation(trench , -trench_x/2 , 0 , 0 )

    au_trench = geompy.MakeBoxDXDYDZ(trench_x+20 , trench_y + data.bfm , trench_z+20)
    au_trench = geompy.MakeTranslation(au_trench , -(trench_x+20)/2 , 0 , 0 )
    
    ridge_trench = geompy.MakeBoxDXDYDZ(22 , base_y, e_stop)
    ridge_trench = geompy.MakeTranslation(ridge_trench , -11 , 0 , 0)

    
    if extra_width:
        extra_width = 1
    else:
        extra_width = 0
    

    base = geompy.MakeBoxDXDYDZ(base_x * (n_active + extra_width) , base_y + data.bfm, base_z)
    base = geompy.MakeTranslation(base , (-base_x / 2) * extra_width , 0 , 0 )

    base_0 = geompy.MakeBoxDXDYDZ(base_x * (n_active + extra_width) , base_y + data.bfm, base_z - 30)
    base_0 = geompy.MakeTranslation(base_0, (-base_x / 2) * extra_width , 0 , 0 )    
        

    for i in range( 0 , n_active): #insert all the active regions into the chip
        x_pos = base_x * ( i + 0.5) 
        ridge_cut = geompy.MakeTranslation(ridge , x_pos , 0 , z_pos)
        ridge_trench_cut = geompy.MakeTranslation(ridge_trench, x_pos, 0, base_z - e_stop)
        base = geompy.MakeCut(base , ridge_cut , True)
        base = geompy.MakeCut(base , ridge_trench_cut , True)
        partition_array = np.append(partition_array, ridge_cut)
    
    #write the body details of the case.sif file for the ridges
    write_ridge_bodies(data, pro_name)
    
    count = n_active * data.n_layers + 1 

    for i in range ( 0 , n_active+1): #create trenches halfway between each active region
        x_pos = base_x * i
        trench_cut = geompy.MakeTranslation(trench , x_pos , 0 , base_z - trench_z)
        base = geompy.MakeCut(base , trench_cut , True)

    if data.bfm !=0: #create a back facet trench if there is a back facet monitor present
        bfm_trench_z = trench_z
        back_trench = geompy.MakeBoxDXDYDZ(base_x * (n_active + 0) , 20 , bfm_trench_z+10)
        back_facet_trench = geompy.MakeTranslation(back_trench, 0, base_y , base_z -  bfm_trench_z )
        base = geompy.MakeCut(base , back_facet_trench , True)
        base_0 = geompy.MakeCut(base_0, back_facet_trench, True)

        back_trench = geompy.MakeBoxDXDYDZ(base_x * (n_active + 1) , 20 , trench_z+10)
        back_facet_trench = geompy.MakeTranslation(back_trench, -base_x/2 , base_y + data.bfm -20, base_z - trench_z )
        base = geompy.MakeCut(base , back_facet_trench , True)
        base_0 = geompy.MakeCut(base_0, back_facet_trench, True)

    
    base = geompy.MakeCut(base, base_0) #separate upper and main part of the chip as a viscous layer

    if data.pside_down !=0: 
        partition_array, count = create_pdown_layer(data, geompy, base, base_0, e_stop, pro_name, count, partition_array)

    count = write_partition(count, data.device_mat, pro_name, name = 'Base buffer')
    partition_array = np.append(partition_array, base)

    if data.current_model == 0:
        count = write_partition(count, data.device_mat, pro_name, name = 'Base lower')
        partition_array = np.append(partition_array, base_0) 

    #add submount or external_heatsink   

    if data.ext_sink_mat !=0:
        partition_array, count = create_submount(data, geompy, partition_array, count, pro_name)

    final_partition = geompy.MakePartition(partition_array.tolist(),  [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

    return(final_partition , ridge, count)   



def create_submesh_group(geompy, fin_partition, explode_partiton, inds):
    #primes the necessary bodies into an 'auto_group' for
    #precise obect parameters/ submesh
    objs = [] 
    for i in inds:
        objs.append(explode_partiton[i])
    sub_mesh_auto_group = geompy.CreateGroup(fin_partition,geompy.ShapeType["SOLID"] )
    geompy.UnionList(sub_mesh_auto_group, objs)
    return(sub_mesh_auto_group)

    

def NETGEN_submesh(sc_data, mesh_11, sm_autgroup, sink = False, local_L = 0.3):
    ##define algo for 1D,2D
    if sink:
        NETGEN_1D_2D_3D_1 = mesh_11.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=sm_autgroup)
        sm_object = NETGEN_1D_2D_3D_1.GetSubMesh()
        NETGEN_3D_Parameters_2 = NETGEN_1D_2D_3D_1.Parameters()
        NETGEN_3D_Parameters_2.SetMaxSize( sc_data.bdy_mesh )
        NETGEN_3D_Parameters_2.SetMinSize( 10)
        NETGEN_3D_Parameters_2.SetSecondOrder( 0 )
        NETGEN_3D_Parameters_2.SetOptimize( 1 )
        NETGEN_3D_Parameters_2.SetFineness( 3)
        NETGEN_3D_Parameters_2.SetChordalError( -1 )
        NETGEN_3D_Parameters_2.SetChordalErrorEnabled( 0 )
        NETGEN_3D_Parameters_2.SetUseSurfaceCurvature( 1 )
        NETGEN_3D_Parameters_2.SetFuseEdges( 1 )
        NETGEN_3D_Parameters_2.SetQuadAllowed( 0 )
        NETGEN_3D_Parameters_2.SetCheckChartBoundary( 176 )

    else:
        NETGEN_1D_2D_3D_2 = mesh_11.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=sm_autgroup)
        sm_object = NETGEN_1D_2D_3D_2.GetSubMesh()
        NETGEN_3D_Simple_Parameters_2 = NETGEN_1D_2D_3D_2.Parameters(smeshBuilder.SIMPLE)
        NETGEN_3D_Simple_Parameters_2.SetLocalLength(sc_data.r_mesh)
        NETGEN_3D_Simple_Parameters_2.LengthFromEdges()
        NETGEN_3D_Simple_Parameters_2.LengthFromFaces()

    return(mesh_11, sm_object)

def NETGEN_create_mesh(sc_data, fin_partition, smesh):

    bdy_dim = sc_data.bdy_mesh ; r_dim = sc_data.r_mesh
    #ridge mesh properties 
    Mesh_1 = smesh.Mesh(fin_partition,'Mesh_1')
    NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
    NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
    NETGEN_3D_Parameters_1.SetMaxSize( 10 )
    NETGEN_3D_Parameters_1.SetMinSize( 0.2)
    NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
    NETGEN_3D_Parameters_1.SetOptimize( 1 )
    NETGEN_3D_Parameters_1.SetFineness( 4 )
    NETGEN_3D_Parameters_1.SetChordalError( -1 )
    NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
    NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
    NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
    NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
    NETGEN_3D_Parameters_1.SetCheckChartBoundary( 176 )
    return(Mesh_1)



def add_to_study(geompy, O_array , mult_ridg, part, part_expl):
    O , OX, OY,OZ = O_array
    geompy.addToStudy( O, 'O' )
    geompy.addToStudy( OX, 'OX' )
    geompy.addToStudy( OY, 'OY' )
    geompy.addToStudy( OZ, 'OZ' )
    geompy.addToStudy(mult_ridg, 'multi ridge' )
    geompy.addToStudy(part , 'final partition')

    for i in range( 0 , len(part_expl)):
        geompy.addToStudyInFather( part, part_expl[i], f'Solid_{i}' )
    return()


def find_face_of_god(smesh , group_array, data, scc  = False):
    
    #find the index of the face for the sink 
    #this is required as there is not consistency for the numbering of each face with respect to creation order
    temp_z = 1
    temp_pz = 1

    if data.pside_down == 0:
        for i  in range(0 , len(group_array)):
                BBox = smesh.GetBoundingBox(group_array[i])
                if BBox.minZ == BBox.maxZ: 

                    if BBox.minZ < temp_z:
                        temp_z = BBox.minZ
                        f_o_g = i + 1
                    if BBox.maxZ > temp_pz and scc:
                        temp_pz = BBox.MaxZ
        return(int(f_o_g))
    
    elif data.pside_down ==1:
        for i  in range(0 , len(group_array)):
                BBox = smesh.GetBoundingBox(group_array[i])
                if BBox.minZ == BBox.maxZ: 
                    if BBox.maxZ > temp_z:
                        temp_z = BBox.maxZ
                        f_o_g = i + 1
                    if BBox.maxZ > temp_pz and scc:
                        temp_pz = BBox.MaxZ
        return(int(f_o_g))
    
def find_wirebonds_of_god(smesh , group_array, data):
    z_thresh = data.device_dim[2] + data.au_cap + data.insul_z+1
    SCC_faces = np.array([])
    for i in range(len(group_array)):
        BBox = smesh.GetBoundingBox(group_array[i])
        if BBox.minZ == BBox.maxZ and BBox.minZ > z_thresh:
            SCC_faces = np.append(SCC_faces, i+1)
    SCC_faces.astype(int)
    return(SCC_faces)

def create_solid_array(Mesh_1 , partition_exploded):

    solid_array = []
    for i in range( 0 , len(partition_exploded)):
        solid_array = np.append(solid_array , Mesh_1.GroupOnGeom(partition_exploded[i],'Solid_0',SMESH.VOLUME)) 
    solid_array =  Mesh_1.GetGroups()
    group_array = Mesh_1.GetMesh().FaceGroupsSeparatedByEdges( 30, 0, 0 ) #seperated by edges for the mesh

    return(solid_array , group_array)

def new_mesh_ext_sink(data, arg0 ): # (ridge mesh , body mesh)

    print(f'SALOME SCRIPT COMMENCING:')

    geompy = geomBuilder.New()
    O = geompy.MakeVertex(0, 0, 0)
    OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
    OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
    OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

    #multi_ridge = create_ridges (data , geompy)
    adams_partition, multi_ridge, count = create_device( data ,  geompy , back_facet_trench=False, middle_y=False, extra_width= True, pro_name=arg0)
    
    count = count -1

    partition_exploded = geompy.ExtractShapes(adams_partition, geompy.ShapeType["SOLID"], False)#exploded the object into a large array
    


    #################-----------------------------------------------------------MESHING-------------------------------------------------
    smesh = smeshBuilder.New()



    if data.pside_down ==0:
        ridge_inds = np.arange(0,data.n_ridges * data.n_layers,1).tolist()
        ridge_smag = create_submesh_group(geompy, adams_partition, partition_exploded, ridge_inds)

        if data.insul_mat !=0 and data.au_cap_mat !=0:
            au_inds , insul_inds = get_wirebond_index(data) 
            au_smag = create_submesh_group(geompy, adams_partition, partition_exploded, au_inds)
            ins_smag = create_submesh_group(geompy, adams_partition, partition_exploded, insul_inds)

            if data.ext_sink_mat !=0:
                sink_smag = create_submesh_group(geompy, adams_partition, partition_exploded, [-2,-1])
            else:
                sink_smag = create_submesh_group(geompy, adams_partition, partition_exploded, [-1])
            #start creating the mesh and submeshes
            Mesh_1 = NETGEN_create_mesh(data, adams_partition, smesh)
            Mesh_1, au_smObj = NETGEN_submesh(data, Mesh_1, au_smag, sink = False )
            Mesh_1, ins_smObj =  NETGEN_submesh(data, Mesh_1, ins_smag, sink = False )
            Mesh_1, sink_smObj =  NETGEN_submesh(data, Mesh_1, sink_smag, sink = True )
            Mesh_1, ridge_smObj = NETGEN_submesh(data, Mesh_1, ridge_smag, sink= False)
            isDone = Mesh_1.SetMeshOrder( [ [ au_smObj, ins_smObj, sink_smObj, ridge_smObj] ])
            #isDone = Mesh_1.Compute()
        else:

            sink_smag = create_submesh_group(geompy, adams_partition, partition_exploded, [-3,-1])
            Mesh_1 = NETGEN_create_mesh(data, adams_partition, smesh)
            Mesh_1, ridge_smObj = NETGEN_submesh(data, Mesh_1, ridge_smag, sink= False)
            Mesh_1, sink_smObj =  NETGEN_submesh(data, Mesh_1, sink_smag, sink = True )
            
            isDone = Mesh_1.SetMeshOrder( [ [ sink_smObj, ridge_smObj] ])
            #isDone = Mesh_1.Compute()


    elif data.pside_down ==1:
        ridge_inds = np.arange(0,count - 3,1).tolist()
        ridge_smag = create_submesh_group(geompy, adams_partition, partition_exploded, ridge_inds)
        sink_smag = create_submesh_group(geompy, adams_partition, partition_exploded, [-3,-1])
        Mesh_1 = NETGEN_create_mesh(data, adams_partition, smesh)
        Mesh_1, sink_smObj =  NETGEN_submesh(data, Mesh_1, sink_smag, sink = True )
        Mesh_1, ridge_smObj = NETGEN_submesh(data, Mesh_1, ridge_smag, sink= False)
        isDone = Mesh_1.SetMeshOrder( [ [ sink_smObj, ridge_smObj] ])


    #isDone = Mesh_1.Compute()


    add_to_study(geompy, [O,OX,OY,OZ] , multi_ridge, adams_partition, partition_exploded)

    
    netgen = True

    solid_array , group_array = create_solid_array(Mesh_1 , partition_exploded)

    fog_fog = find_face_of_god(smesh , group_array, data) 
    if data.current_model ==1:
        wb_fog = find_wirebonds_of_god(smesh, group_array, data)


    try:
      Mesh_1.ExportUNV( f'C:/ElmerFEM/ElmerFEM/bin/{arg0}/temp_save.unv', 0 )
      pass
    except:
      print('ExportUNV() failed. Invalid file name?')

    print('#####################################################')
    print(pd.__file__)
    print('#####################################################')

    if data.current_model ==0:
        return(fog_fog)
    elif data.current_model ==1:
        return(fog_fog, wb_fog)
    


arg1 = 'y' #project name - can also take in more arguments if necessary

V = 0 #input variable
sweeping_V = 0


pandas_data = pd.read_csv('C:/Projects/Perry_run/input_csv.csv').to_numpy()
device = MySemiconductor(pandas_data) #put data into MySemiconductor class
device.device_arb_parameter = -1 #placeholder to have the device ridge slightly not off the edge of the whole thing

device.pside_down = 1
device.current_model = 0
device.insul_z = 0.3
device.au_cap = 2
device.ext_sink_dim = [2000, 1500, 300]


new_mesh_ext_sink(device, arg1)

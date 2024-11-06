import numpy as np
import os
import sys

import_variable = 101

instance_to_transfer = 14 #:)




class MySemiconductor:
    def __init__(self, input_dat):
        #device
        self.n_ridges = int(input_dat[1, 0 ]) #number of ridges
        self.r_onoff = input_dat[1:1+int(input_dat[1, 0 ]),1].astype(int)
        self.device_dim = input_dat[1:4,2].astype(float) #physical dimensions of the device (per ridge)
        self.device_mat = int(input_dat[1,3])
        self.n_layers = int(input_dat[1,4])
        #ridges
        self.r_heights =input_dat[1:1+int(input_dat[1,4]),5].astype(float)
        self.r_widths = input_dat[1:1+int(input_dat[1,4]),6].astype(float)
        self.r_lengths = input_dat[1:1+int(input_dat[1,4]),7].astype(float)
        #below are details only required for elmer
        self.r_materials = input_dat[1:1+int(input_dat[1,4]),8].astype(int)
        self.r_heat_power = input_dat[1:1+int(input_dat[1,4]),9].astype(float)
        self.T_sink = float(input_dat[1,10])
        self.bdy_mesh = float(input_dat[1,11])
        self.r_mesh = float(input_dat[1,12])
        self.z_ridge = float(input_dat[1,13])
        self.ext_sink_dim = input_dat[1:4, 14].astype(float)
        self.ext_sink_mat = int(input_dat[1,15])
        #front and back facet details
        self.au_cap = float(input_dat[1,16])
        self.au_cap_mat = float(input_dat[1,17])
        self.bfm = float(input_dat[1,18])
        self.thermistor_dim = input_dat[1:4,19].astype(float)
        self.thermistor_mat = int(input_dat[1,20])

        self.device_arb_parameter = float(0) #this is an arbitrary parameter, used as a placeholder to sweep something directly from the code
        self.trench_dim = input_dat[1:4 , 21].astype(float)
        self.cartridge_mat = int(input_dat[1,22])
        self.cartridge_dim = input_dat[1:4,23].astype(float)
        self.n_chips = int(input_dat[1, 24])

        self.insul_mat = int(input_dat[1,1])
        self.insul_z = float(input_dat[1,1])

def write_header(proj_name,file):
    header = f'Header\n  CHECK KEYWORDS Warn\n  Mesh DB "{proj_name}" "dummy"\n  Include Path ""\n  Results Directory ""\nEnd\n\n'
    file.write(header)
    return()

def write_simconst(file):
    sim = 'Simulation\n  Max Output Level = 5\n  Coordinate System = Cartesian\n  Coordinate Mapping(3) = 1 2 3\n  Simulation Type = Steady state\n  Steady State Max Iterations = 1\n  Output Intervals(1) = 1\n  Coordinate Scaling = 1e-6\n  Solver Input File = case.sif\n  Post File = case.vtu\nEnd\n\nConstants\n  Gravity(4) = 0 -1 0 9.82\n  Stefan Boltzmann = 5.670374419e-08\n  Permittivity of Vacuum = 8.85418781e-12\n  Permeability of Vacuum = 1.25663706e-6\n  Boltzmann Constant = 1.380649e-23\n  Unit Charge = 1.6021766e-19\nEnd\n'
    file.write(sim)
    return()

def write_ind_body(ind , material , body_force): # write an individual body force term
    #if body force is zero then remove body force from sif text
    #material and body_force are integers which indicate the index of the material and forces later described in the .sif file. 
    if body_force ==0:
        bdy_string00 = f'Body {ind} \n  Target Bodies(1) = {ind}  \n  Name = "Body Property {ind}" \n  Equation = 1 \n  Material = {material} \n  Initial condition = 1 \nEnd \n'
    else:
        bdy_string00 = f'Body {ind} \n  Target Bodies(1) = {ind}  \n  Name = "Body Property {ind}" \n  Equation = 1 \n  Material = {material} \n  Body Force = {body_force} \n  Initial condition = 1 \nEnd \n'
    return (bdy_string00)

def write_bodies( device,file): 
    body_string = '\n'

    n_r = device.n_ridges #number of ridges per chip
    n_l = device.n_layers #number of layers per ridge
    if device.thermistor_mat !=0:
        t_count = 1
    else:
        t_count = 0

    if device.cartridge_mat!=0: #mulitple chips - therefore cartridge and chuck
        n_c = device.n_chips #number of chips
        count = 1 #counter for indexing the body force (heating power)
        
        for c in range(0, n_c): #chip index
            for r in range(0, n_r): #ridge index
                count = 1
                for l in range(1,n_l+1): #layer index
                    num = c*(3 + t_count + n_r*n_l) + r*(n_l) + l
                    if device.r_heat_power[l-1] *device.r_onoff[r] !=0: #this means there IS heating in this layer
                        body_string = body_string + write_ind_body(num, device.r_materials[l-1] ,count) + '\n'
                        count += 1
                    else:
                        body_string = body_string + write_ind_body(num, device.r_materials[l-1] , 0) + '\n'
            

            num = (c+1)*(n_r*n_l + 3 + t_count) - 2 - t_count

            body_string = body_string + write_ind_body(num, device.device_mat , 0) + '\n' #chip base
            body_string = body_string + write_ind_body(num+1, device.ext_sink_mat , 0) + '\n'#chip submount
            body_string = body_string + write_ind_body(num+2, 9 , 0) + '\n'#submount thermal paste

            if t_count ==1:
                body_string = body_string + write_ind_body(num + 2 + t_count, device.thermistor_mat , 0) + '\n'#thermistor on submount

        body_string = body_string + write_ind_body(num+3 + t_count , device.cartridge_mat , 0 ) + '\n' #cartridge
        body_string = body_string + write_ind_body(num+4 + t_count , device.cartridge_mat , 0 ) + '\n' #chuck sidewall
        body_string = body_string + write_ind_body(num+5 + t_count , device.cartridge_mat , 0 ) + '\n' #chuck sidewall
        body_string = body_string + write_ind_body(num+6 + t_count , device.cartridge_mat , 0 ) + '\n' #chuck base

    else: #single chip - therefore no cartridge and chuck
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
        

        if device.au_cap_mat != 0 :
            for i in range(int(3*n_r)):
                body_string = body_string + write_ind_body(num + i, int(device.au_cap_mat), 0)  
            num += int(3*(n_r))

        body_string = body_string + write_ind_body(num , device.device_mat , 0) + f'\n' #chip base viscouse upper
        body_string = body_string + write_ind_body(num +1 , device.device_mat , 0) + f'\n' #chip base
        num +=2
        if device.ext_sink_mat !=0:
            body_string = body_string + write_ind_body(num, device.ext_sink_mat , 0) + '\n'#chip submount
            #body_string = body_string + write_ind_body(num+1, 9 , 0) + '\n'#chip submount
            num +=2
        if t_count ==1:
                body_string = body_string + write_ind_body(num + 1 + t_count, device.thermistor_mat , 0) + '\n'#thermistor on submount
        
    file.write(f'{body_string}')
    return()


def write_solver(file): #write the solver section - requires not specific input. Parameters here can be altered for more advanced simulations
    #solver = ' \n Solver 1 \n  Equation = Navier-Stokes\n  Procedure = "FlowSolve" "FlowSolver"\n  Variable = Flow Solution[Velocity:3 Pressure:1]\n  Exec Solver = Always \n  Stabilize = True \n  Optimize Bandwidth = True \n  Steady State Convergence Tolerance = 1.0e-5 \n  Nonlinear System Convergence Tolerance = 1.0e-7 \n  Nonlinear System Max Iterations = 20 \n  Nonlinear System Newton After Iterations = 3 \n  Nonlinear System Newton After Tolerance = 1.0e-3 \n  Nonlinear System Relaxation Factor = 1 \n  Linear System Solver = Iterative \n  Linear System Iterative Method = BiCGStab \n  Linear System Max Iterations = 500 \n  Linear System Convergence Tolerance = 1.0e-10 \n  BiCGstabl polynomial degree = 2 \n  Linear System Preconditioning = ILU0 \n  Linear System ILUT Tolerance = 1.0e-3 \n  Linear System Abort Not Converged = False \n  Linear System Residual Output = 10 \n  Linear System Precondition Recompute = 1 \n End \n Solver 2\n  Equation = Heat Equation\n  Variable = Temperature\n  Procedure = "HeatSolve" "HeatSolver"\n  Exec Solver = Always\n  Stabilize = True\n  Optimize Bandwidth = True\n  Steady State Convergence Tolerance = 1.0e-5\n  Nonlinear System Convergence Tolerance = 1.0e-7\n  Nonlinear System Max Iterations = 20\n  Nonlinear System Newton After Iterations = 3\n  Nonlinear System Newton After Tolerance = 1.0e-3\n  Nonlinear System Relaxation Factor = 1\n  Linear System Solver = Iterative\n  Linear System Iterative Method = BiCGStab\n  Linear System Max Iterations = 500\n  Linear System Convergence Tolerance = 1.0e-10\n  BiCGstabl polynomial degree = 2\n  Linear System Preconditioning = ILU0\n  Linear System ILUT Tolerance = 1.0e-3\n  Linear System Abort Not Converged = False\n  Linear System Residual Output = 10\n  Linear System Precondition Recompute = 1\nEnd'
    solver1 = open('C:/Projects/Perry_run/case_text/case_solver.txt' , 'r')
    solver = solver1.read()
    file.write(solver)
    solver1.close()
    return()

def write_equation(file): #write the equation section - simply activates the heat solver within the multiphysics simulations
    eqn = '\n\nEquation 1\n  Name = "Heat equation"\n  Active Solvers(1) = 1\nEnd\n\n'
    file.write(eqn)
    return()

def write_materials(file): #write materials and their properties - this can be edited quite easily and perhaps could be made to include new materials as an input?
    #mats = 'Material 1\n  Name = "XX Indium Phosphide"\n  Heat Conductivity = 68\n  Heat Capacity = 310\n  Density = 4800\n  Electric Conductivity = 10e6\n  Porosity Model = Always saturated\nEnd\n\nMaterial 2\n  Name = "XX Gold"\n  Porosity Model = Always saturated\n  Heat Capacity = 129\n  Heat Conductivity = 318\n  Density = 19300\n  Electric Conductivity = 4.1e7\nEnd\n\nMaterial 3\n  Name = "XX Silicon Nitride"\n  Porosity Model = Always saturated\n  Density = 3230\n  Heat Capacity = 880\n  Electric Conductivity = 10e-14\n  Heat Conductivity = 17.5\nEnd\n\nMaterial 4\n  Name = "XX Titanium"\n  Heat Capacity = 522.4\n  Heat Conductivity = 11.4\n  Electric Conductivity = 2.38e6\n  Porosity Model = Always saturated\n  Density = 4510\nEnd\n\nMaterial 5\n  Name = "XX Platinum"\n  Heat Capacity = 134\n  Heat Conductivity = 69.1\n  Electric Conductivity = 9.43e6\n  Density = 21450\nEnd\n\nMaterial 6\n  Name = "XX AlN"\n  Heat Capacity = 740\n  Density = 2920\n  Porosity Model = Always saturated\n  Heat Conductivity = 237\nEnd\n\nMaterial 7\n  Name = "XX Ta2O5"\n  Heat Capacity = 300\n  Density = 8200\n  Porosity Model = Always saturated\n  Heat Conductivity = 7\nEnd\n\n'
    mats_1 = open("C:/Projects/Perry_run/case_text/case_materials.txt" , 'r')
    mats = mats_1.read()
    file.write(mats)
    mats_1.close()
    return() 

def write_body_forces(device, file):
    heat_power = device.r_heat_power
    body_on = device.r_onoff
    if device.cartridge_mat!=0: #multiple chips
        multiplier = np.sum(body_on) * device.n_chips
    else: #single chip
        multiplier = np.sum(body_on) 
    count = 1
    for i in heat_power:
        if i != 0:
            bdy_force = f'Body Force {count}\n  Name = "Body Force {count}"\n  Integral Heat Source = {i * multiplier}\n  Heat Source = 1\nEnd\n\n'
            file.write(bdy_force)
            count += 1
        else:
            pass
    return() 

def write_initial_conds(file ,T_init = 60): #initial temp of all objects, doesnt matter too much for steady state simulations 
    init_c = f'Initial Condition 1\n  Name = "Initial Temperature"\n  Temperature = {T_init}\nEnd'
    file.write(init_c)
    return()

def write_boundary_conds(boundary , file ,T_sink = 80): #heat sink temperature - though is not necessarily that complicated to change the input into an arry
    #here we need the boundary number of the lowest boundary (which was far too much effort to determing from the stupid Salome simulations)
    bound_cond = f'Boundary Condition 1\n  Target Boundaries(1) = {boundary}\n  Name = "Heat Sink"\n  Temperature = {T_sink}\nEnd'
    file.write(bound_cond)
    return()

#write elmer solver start info - deal with this seperately
def write_ESSI(proj_name):#ELMERSOLVER_STARTINFO write
    directory = 'C:/ElmerFEM/ElmerFEM/bin/'
    e = open(directory + 'ELMERSOLVER_STARTINFO', "w")
    e.write(f'{proj_name}\case.sif \n1')
    e.close()
    return()

def global_write(project_name , device):
    my_file = open(f'C:/ElmerFEM/ElmerFEM/bin/{project_name}/case.sif' , 'w')
    write_header(project_name, my_file) #header needs project name
    my_file.close() ; my_file = open(f'C:/ElmerFEM/ElmerFEM/bin/{project_name}/case.sif' , 'a')
    write_simconst(my_file)
    write_bodies(device , my_file)# requires the number of ridges - layers - materials
    write_solver(my_file)
    write_equation(my_file)
    write_materials(my_file)
    write_body_forces(device,my_file) #needs heat power for each ridge
    write_initial_conds(my_file , T_init = device.T_sink+10)#can take in the a different initial temp if that is of interest
    #write_boundary_conds -DO NOT EXECUTE HERE, needs Salome input so execute direclty from Salome script
    write_ESSI(project_name)
    my_file.close()
    return()
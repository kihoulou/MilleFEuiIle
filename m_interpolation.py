from dolfin import *
import numpy
import math
import os 
from random import uniform
import sys
from m_parameters import *
from m_material_properties import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

def Max(a, b): return (a+b+abs(a-b))/2.0

def Min(a, b): return (a+b-abs(a-b))/2.0

def tracer_counting(mesh, tracers_in_cells, tracers):
    n1 = 0.0
    n2 = 0.0

    for j in range(mesh.num_cells()):
        for i in range(0,len(tracers_in_cells[j])):
            tracer_no = tracers_in_cells[j][i]
            if (tracers[tracer_no][11][0] == 1.0): # denser
                n1 += 1.0
            if (tracers[tracer_no][11][1] == 1.0): # lighter
                n2 += 1.0

    return n1, n2

 

def composition_interpolation(mesh, tracers_in_cells, tracers, column, function):    
# def composition_interpolation(mesh, tracers_in_cells, tracers, column, function, n_tracers_orig, n_tracers_current):    

    ranks = []
    for j in range(mesh.num_cells()):
        weight_cell = 0.0
        value_cell = 0.0

        if (len(tracers_in_cells[j]) > 0):
            for i in range(0, len(tracers_in_cells[j])):
                tracer_no = tracers_in_cells[j][i]
            
                # if (weight_tracers_by_ratio == True):
                #     weight_tracer = n_tracers_orig[column]/n_tracers_current[column]
                # else:
                weight_tracer = 1.0
                
                value_cell += weight_tracer*tracers[tracer_no][10][column]
                weight_cell += weight_tracer

            ranks.append(value_cell/weight_cell)
        else:
            if (column == empty_cells_composition):
                ranks.append(1.0)
            else:
                ranks.append(0.0)

    ranks = numpy.array(ranks)                                                            
    function.vector().set_local(ranks)


def melt_interpolation(mesh, tracers_in_cells, tracers, column, function):    
    ranks = []
    for j in range(mesh.num_cells()):
        weight_cell = 0.0
        value_cell = 0.0

        if (len(tracers_in_cells[j])==0):
            ranks.append(0.0)
        
        else:
            for i in range(0,len(tracers_in_cells[j])):
                tracer_no = tracers_in_cells[j][i]
                weight_tracer = 1.0
                
                value_cell += weight_tracer*tracers[tracer_no][11][column]
                weight_cell += weight_tracer

            ranks.append(value_cell/weight_cell)

    ranks = numpy.array(ranks)                                                            
    function.vector().set_local(ranks)

def scalar_interpolation(mesh, tracers_in_cells, tracers, column, avg_type, function):
    ranks = []
    for j in range(mesh.num_cells()):
        weight_cell = 0.0
        value_cell = 0.0
        centroid = Cell(mesh, j).midpoint()

        if (len(tracers_in_cells[j])==0):
            # print("Interpolation: Zero tracers in cell", j, "len PIC",len(tracers_in_cells[j]))
            # print("Cell coordinates:", centroid.x(), centroid.y())
            # print("Prescribing zero...")
            ranks.append(0.0)

            # print("Exiting...")
            # MPI.comm_world.Abort()

        if (avg_type=="ARITM" and len(tracers_in_cells[j])!=0):
            for i in range(0,len(tracers_in_cells[j])):
                tracer_no = tracers_in_cells[j][i]
                value_cell += tracers[tracer_no][column]
                weight_cell += 1.0 

            ranks.append(value_cell/weight_cell)

        if (avg_type=="HARM" and len(tracers_in_cells[j])!=0 ):
            for i in range(0,len(tracers_in_cells[j])):
                tracer_no = tracers_in_cells[j][i]
                value_cell += ln(tracers[tracer_no][column])/ln(10.0)
                weight_cell += 1.0 

            ranks.append(10.0**(value_cell/weight_cell))
        
    ranks = numpy.array(ranks)                                                            
    function.vector().set_local(ranks)

def tracer_count_interpolation(mesh, tracers_in_cells, function):
    ranks = []
    for j in range(mesh.num_cells()):
        ranks.append(len(tracers_in_cells[j]))

    ranks = numpy.array(ranks)                                                            
    function.vector().set_local(ranks)

def stress_interpolation(mesh, tracers_in_cells, tracers, stress_tensor):
    ranks = []

    for j in range(mesh.num_cells()):
        weight_cell = 0.0

        if (len(tracers_in_cells[j]) == 0):
            ranks.append(0.0)
            ranks.append(0.0)
            ranks.append(0.0)
            ranks.append(0.0)
            
        # --- Linear averaging of the stress ---
        # t_xx = 0.0
        # t_xz = 0.0

        # if (len(tracers_in_cells[j])!=0):
        #     for i in range(0,len(tracers_in_cells[j])):
        #         tracer_no = tracers_in_cells[j][i]

        #         t_xx += tracers[tracer_no][4]
        #         t_xz += tracers[tracer_no][5]
        #         weight_cell += 1.0 

        #     ranks.append(t_xx/weight_cell)
        #     ranks.append(t_xz/weight_cell)
        #     ranks.append(t_xz/weight_cell)
        #     ranks.append(-t_xx/weight_cell)

        # --- Quasi-geometric averaging of the stress ---
        t_xx_pos = 0.0
        t_xx_neg = 0.0

        t_xz_pos = 0.0
        t_xz_neg = 0.0
        
        if (len(tracers_in_cells[j])!=0):
            for i in range(0,len(tracers_in_cells[j])):
                tracer_no = tracers_in_cells[j][i]
                
                if (tracers[tracer_no][4] > 0.0):
                    t_xx_pos += ln(tracers[tracer_no][4])/ln(10.0)
                elif (tracers[tracer_no][4] < 0.0):
                    t_xx_neg += ln(-tracers[tracer_no][4])/ln(10.0)
                
                if (tracers[tracer_no][5] > 0.0):
                    t_xz_pos += ln(tracers[tracer_no][5])/ln(10.0)
                elif (tracers[tracer_no][5] < 0.0):
                    t_xz_neg += ln(-tracers[tracer_no][5])/ln(10.0)

                weight_cell += 1.0 

            ranks.append(0.5*(10**(t_xx_pos/weight_cell) - 10**(t_xx_neg/weight_cell)))
            ranks.append(0.5*(10**(t_xz_pos/weight_cell) - 10**(t_xz_neg/weight_cell)))
            ranks.append(0.5*(10**(t_xz_pos/weight_cell) - 10**(t_xz_neg/weight_cell)))
            ranks.append(-0.5*(10**(t_xx_pos/weight_cell) - 10**(t_xx_neg/weight_cell)))

    ranks = numpy.array(ranks) 
    stress_tensor.vector().set_local(ranks)

def stress_update(mesh, tracers_in_cells, tracers, visc, strain_rate_tensor, z_function):
    for j in range(mesh.num_cells()):
        for i in range(0, len(tracers_in_cells[j])):
            tracer_no = tracers_in_cells[j][i]

            px = tracers[tracer_no][0]
            pz = tracers[tracer_no][1]

            z_tracer = z_function(Point(px, pz))
            visc_tracer = visc(Point(px, pz))

            e_xx = strain_rate_tensor(Point(px, pz))[0]
            e_xz = strain_rate_tensor(Point(px, pz))[1]
            e_zz = strain_rate_tensor(Point(px, pz))[3]

            tracers[tracer_no][4]  += z_tracer*(2*visc_tracer*e_xx - tracers[tracer_no][4])
            tracers[tracer_no][5]  += z_tracer*(2*visc_tracer*e_xz - tracers[tracer_no][5])
 
def stress_reduction(mesh, tracers_in_cells, tracers, yield_function, yield_stress, stress_invariant):
    for j in range(mesh.num_cells()):
        for i in range(0,len(tracers_in_cells[j])):
            tracer_no = tracers_in_cells[j][i]

            px = tracers[tracer_no][0]
            py = tracers[tracer_no][1]

            yield_function_tracer = yield_function(Point(px,py))
            
            if (yield_function_tracer >= -1e-3):
                yield_stress_tracer = yield_stress(Point(px,py))
                stress_invariant_tracer = stress_invariant(Point(px,py))

                tracers[tracer_no][4]  = tracers[tracer_no][4]*yield_stress_tracer/stress_invariant_tracer
                tracers[tracer_no][5]  = tracers[tracer_no][5]*yield_stress_tracer/stress_invariant_tracer

def stress_rotation(mesh, tracers_in_cells, tracers, vorticity_tensor, dt):
    for j in range(mesh.num_cells()):
        for i in range(0,len(tracers_in_cells[j])):
            tracer_no = tracers_in_cells[j][i]

            px = tracers[tracer_no][0]
            pz = tracers[tracer_no][1]

            t_xx = tracers[tracer_no][4]
            t_xz = tracers[tracer_no][5]

            W_xz = vorticity_tensor(Point(px,pz))[1]

            # t_old = t - dt*(t.W - W.t) 
            tracers[tracer_no][4] += - float(dt)*(2*t_xz*W_xz)
            tracers[tracer_no][5] += + float(dt)*(2*t_xx*W_xz)

from dolfin import *
from m_parameters import *
from m_material_properties import *
from m_interpolation import *

class Melting:

    def __init__(self, MeshClass, ElemClass, TracersClass):
        
        self.mesh = MeshClass.mesh
        
        self.Temp = ElemClass.Temp
        self.composition = ElemClass.composition
        self.xm = ElemClass.xm
        self.xm_k = ElemClass.xm_k

        if (internal_melting == True):
            self.tracers_in_cells = TracersClass.tracers_in_cells
            self.vacancy = TracersClass.vacancy
            self.tracers = TracersClass.tracers
            self.only_melt_tracers = TracersClass.only_melt_tracers

            if (reload_tracers == False):
                    if (self.only_melt_tracers == True):
                        self.add_melt_tracers()

    def add_melt_tracers(self):

        for j in range(self.mesh.num_cells()):
            centroid = Cell(self.mesh, j).midpoint()
        
            r_in = Cell(self.mesh, j).inradius()
            deg2rad = np.pi/180.0

            temp_cell = self.Temp(Point(centroid.x(),centroid.y()))
            xm_cell = self.xm(Point(centroid.x(),centroid.y()))

            if (temp_cell > T_melt - dT_melt_tresh  and len(self.tracers_in_cells[j]) < 10):
                for ii in range(0,5):
                    tracer_angle = 90.0*(ii-1)

                    if (ii==0):
                        # First tracer to the cell center
                        xx = centroid.x()
                        yy = centroid.y()
                    else:
                        # Other tracers along a circle
                        xx = centroid.x() + r_in*0.5*sin(tracer_angle*deg2rad)
                        yy = centroid.y() + r_in*0.5*cos(tracer_angle*deg2rad)

                    if (len(self.vacancy) > 0): # if possible, put them in the vacancies
                        self.tracers[self.vacancy[0]] = [xx, yy, j, rank, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, [xm_cell, 0.0], 1, 0]
                        self.tracers_in_cells[j].append(self.vacancy[0])
                        del self.vacancy[0]

                    else: # if no vacancy is available, put them at the end of the list
                        self.tracers.append([xx, yy, j, rank, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, [xm_cell, 0.0], 1, 0])
                        self.tracers_in_cells[j].append(len(self.tracers)-1) 

    # When the temperature is below melting point and no melt is present, delete these tracers
    def delete_melt_tracers(self):

        MPI.barrier(comm)
        for j in range(self.mesh.num_cells()):
            to_remove = []
            
            for i in range(0,len(self.tracers_in_cells[j])):
                tracer_no = self.tracers_in_cells[j][i]
                xx = self.tracers[tracer_no][0]
                yy = self.tracers[tracer_no][1]

                temp_tracer = self.Temp(Point(xx,yy))
                xm_tracer = self.tracers[tracer_no][11][0]

                if (temp_tracer < T_melt - dT_melt_tresh  and xm_tracer == 0.0):
                    to_remove.append(tracer_no)

                    self.vacancy.append(tracer_no)
                    self.tracers[tracer_no] = []

            to_remove.sort(reverse=True)
            for l in range(len(to_remove)):
                self.tracers_in_cells[j].remove(to_remove[l])


    def compute_melting(self):

        for j in range(self.mesh.num_cells()):
            for i in range(0, len(self.tracers_in_cells[j])):
                tracer_no = self.tracers_in_cells[j][i]
                xx = self.tracers[tracer_no][0]
                yy = self.tracers[tracer_no][1]

                temp_tracer = self.Temp(Point(xx,yy))

                if (temp_tracer >= T_melt):
                    # Ice is melting
                    self.tracers[tracer_no][11][0] += cp(T_melt, self.composition)/Lt*(temp_tracer-T_melt)
                    self.tracers[tracer_no][11][1] = T_melt - temp_tracer

                else:
                    if (self.tracers[tracer_no][11][0] > 0.0):
                        if (abs(cp(T_melt, self.composition)/Lt*(temp_tracer-T_melt)) <= self.tracers[tracer_no][11][0]):
                            # Part of the melt crystallizes
                            self.tracers[tracer_no][11][0] += cp(T_melt, self.composition)/Lt*(temp_tracer-T_melt)
                            self.tracers[tracer_no][11][1] = T_melt - temp_tracer
                        else:
                            # All the melt crystallizes
                            self.tracers[tracer_no][11][1] = self.tracers[tracer_no][11][0]*Lt/cp(T_melt, self.composition) # This one imperatively first !
                            self.tracers[tracer_no][11][0] = 0.0
                            
                    else:
                        # Ice is not melting
                        self.tracers[tracer_no][11][0] = 0.0
                        self.tracers[tracer_no][11][1] = 0.0

    
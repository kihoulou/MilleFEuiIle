from dolfin import *
from m_parameters import *
from m_material_properties import *

# If the temperature approaches the melting point, add tracers to the area


def compute_melting(mesh, tracers_in_cells, tracers, Temp):
    for j in range(mesh.num_cells()):
        for i in range(0,len(tracers_in_cells[j])):
            tracer_no = tracers_in_cells[j][i]
            xx = tracers[tracer_no][0]
            yy = tracers[tracer_no][1]

            temp_tracer = Temp(Point(xx,yy))

            if (temp_tracer >= T_melt):
                # Ice is melting
                tracers[tracer_no][12][0] += cp(T_melt)/Lt*(temp_tracer-T_melt)
                tracers[tracer_no][12][1] = T_melt - temp_tracer

            else:
                if (tracers[tracer_no][12][0] > 0.0):
                    if (abs(cp(T_melt)/Lt*(temp_tracer-T_melt)) <= tracers[tracer_no][12][0]):
                        # Part of the melt crystallizes
                        tracers[tracer_no][12][0] += cp(T_melt)/Lt*(temp_tracer-T_melt)
                        tracers[tracer_no][12][1] = T_melt - temp_tracer
                    else:
                        # All the melt crystallizes
                        tracers[tracer_no][12][1] = tracers[tracer_no][12][0]*Lt/cp(T_melt) # This one imperatively first !
                        tracers[tracer_no][12][0] = 0.0
                        
                else:
                    # Ice is not melting
                    tracers[tracer_no][12][0] = 0.0
                    tracers[tracer_no][12][1] = 0.0
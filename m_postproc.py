from dolfin import *
from m_parameters import *
     
# Root-mean-square velocity
def rms_vel(v):
        return sqrt(assemble(dot(v,v)*dx)/(height*length)) 

# Nusselt number
def nusselt(q_cond_top, q_top):
        return q_top/q_cond_top

# For van Keken (1997)
def entrainment(comp_1):
    phi = Expression("0.0 + 1.0*(x[1] > d_e)", d_e = d_e, degree=1)
    return 1.0/(length*d_b)*assemble(phi*(comp_1)*dx)

# For Fullsack (1995)
def closest_point_interpolation(scheme, points, tracers, vacancy, t):

    for i in range(0,len(points)):
        dist_min= 1.0 # initial value
        j_min = len(tracers)
        for j in range(0,len(tracers)):
            if (j not in vacancy):
                dist = sqrt((points[i][0]-tracers[j][0])**2 + (points[i][1]-tracers[j][1])**2)
                if (dist < dist_min):
                    dist_min = dist
                    j_min = j

        r0 = 0.25
        w0 = 0.30
        
        r = sqrt(points[i][0]**2 + points[i][1]**2)
        w = w0*r/r0*exp(-r/r0)
        
        F =  points[i][0]*cos(w*t) + points[i][1]*sin(w*t)

        file = open("data_"+name+"/tracers/point_"+str(i)+".dat","a")
        file.write((3*"%.6E\t"+"\n")%(float(t), tracers[j_min][5], F))
        file.close()


def melt_height(mesh, tracers, points_in_cells):
    yy = 0.0
    for j in range(mesh.num_cells()):
        for i in range(0,len(points_in_cells[j])):
            tracer_no = points_in_cells[j][i]
            
            if (tracers[tracer_no][12] > 0.0 and tracers[tracer_no][1] > yy):
                yy = tracers[tracer_no][1]

    return yy

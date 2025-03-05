from dolfin import *
import numpy
import os 
from m_parameters import *
import time

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

    
def advect_tracers(moving_tracers, vacancy, tracers_in_cells, name, method, tracers, t, step, dt, v, v_mesh, output_now):
    if (output_now==True and rank==0):
        
        # I do not include composition and melt fraction into the output (zero)
        if (save_tracers == True):
            file = open("data_"+name+"/tracers/step_"+str(step)+".dat","a")
            file.write((2*"%s\t"+"\n")%("#x_position",\
                                        "y_position"))    
            file.close()

        # file.write((11*"%s\t"+"\n")%("#x_position",\
        #                             "y_position",\
        #                             "tau_old_xx",\
        #                             "tau_old_xz",\
        #                             "eps_plast",\
        #                             "ocean_mat",\
        #                             "surface_mat",\
        #                             "orig. depth.",\
        #                             "composition",\
        #                             "melt fraction",\
        #                             "origin"))    
        # file.close()

        # file = open("data_"+name+"/tracers/surface_material/step_"+str(step)+".dat","w")
        # file.close()

        # file = open("data_"+name+"/tracers/ocean_material/step_"+str(step)+".dat","w")
        # file.close()

    MPI.barrier(comm)
    to_remove = []

    # Copy "vacancy" to a separate list, so that it does not change length
    # while in the tracer loop
    vacancy_old = vacancy.copy()
    vac_detected = 0

    for j in range (0,len(tracers)):  

        # We need to exclude vacant positions
        # This one obviously goes through the "vacancy" loop for every tracer --> is quite slow
        # if (j not in vacancy_old):

        # This one is faster: we sort the "vacancy" beforehand and than just skip
        # its components while going through "tracers" loop
            if (len(vacancy_old)>0 and j==vacancy_old[vac_detected]):
                vac_detected = min(vac_detected+1, len(vacancy_old)-1)
                continue 

            px = tracers[j][0]
            py = tracers[j][1]
            v1 = v(Point(px,py))

            try:    
                if (method=="Euler"):
                    tracers[j][0] = px + v1[0]*float(dt)
                    tracers[j][1] = py + v1[1]*float(dt)

                if (method=="RK2"):
                    v2 = v(Point(px+v1[0]*float(dt), py+v1[1]*float(dt)))

                    vx_RK2 = (v1[0] + v2[0])/2.0
                    vy_RK2 = (v1[1] + v2[1])/2.0

                    tracers[j][0] = px + vx_RK2*float(dt)
                    tracers[j][1] = py + vy_RK2*float(dt)

                if (method=="RK4"):
                    # v1 is in the right process and domain by definition
                    # Compute the remaining steps of RK scheme
                    v2 = v(Point(px+v1[0]*float(dt)/2.0,py+v1[1]*float(dt)/2.0))
                    v3 = v(Point(px+v2[0]*float(dt)/2.0,py+v2[1]*float(dt)/2.0))
                    v4 = v(Point(px+v3[0]*float(dt),    py+v3[1]*float(dt)))

                    vx_RK4 = (v1[0]+2.0*v2[0]+2.0*v3[0]+v4[0])/6.0
                    vy_RK4 = (v1[1]+2.0*v2[1]+2.0*v3[1]+v4[1])/6.0

                    # If the final position is in the domain, will it be even after mesh  moves?
                    # Use the mesh velocity in the new position
                    vm = v_mesh(Point(px + vx_RK4*float(dt),py + vy_RK4*float(dt)))
                    
                    # Test for new mesh position
                    v_test = v(Point(px + vx_RK4*float(dt),py + vy_RK4*float(dt)+1.0*vm[1]*float(dt)))
                    v_test = v(Point(px + vx_RK4*float(dt),py + vy_RK4*float(dt)-1.0*vm[1]*float(dt)))

                    # Only then update the tracer position
                    tracers[j][0] = px + vx_RK4*float(dt)
                    tracers[j][1] = py + vy_RK4*float(dt)

                    # !!! The mesh move needs to be checked before updating the position.
                    # If the position is updated, but wouldn't be in the rank after mesh moved,
                    # it would go to the "except part" and would be advected two times !!!

                # Find out the rank at the final (!!) position of the tracer, not earlier!
                v_test = v(Point(tracers[j][0], tracers[j][1]))
                tracers[j][4] = rank

                # Output surface material tracers in every step
                # if (tracers[j][9] > 0.0):
                #     file = open("data_"+name+"/tracers/surface_material/step_"+str(step)+".dat","a")
                #     file.write((2*"%.5E\t"+"%d"+"\n")%(tracers[j][0],tracers[j][1],tracers[j][13]))    
                #     file.close()

                # if (tracers[j][8] > 0.0):
                #     file = open("data_"+name+"/tracers/ocean_material/step_"+str(step)+".dat","a")
                #     file.write((2*"%.5E\t"+"\n")%(tracers[j][0],tracers[j][1]))    
                #     file.close()

                # Output all tracers on dedicated steps
                if (save_tracers == True and output_now == True):


                    file = open("data_"+name+"/tracers/step_"+str(step)+".dat","a")
                    file.write((2*"%.7E\t")%(tracers[j][0], tracers[j][1]))    
        
                    if ("rank" in Tracers_Output):    
                        file.write(("%d\t")%(rank))
                    
                    if ("sigma_xx" in Tracers_Output):    
                        file.write(("%.7E\t")%(tracers[j][5]))

                    if ("sigma_xz" in Tracers_Output):    
                        file.write(("%.7E\t")%(tracers[j][6]))    
                    
                    if ("plastic strain" in Tracers_Output):    
                        file.write(("%.7E\t")%(tracers[j][6]))    

                    

                    # if (tracers[j][11][0] == 1.0): # Output only of the bottom material
                    file = open("data_"+name+"/tracers/step_"+str(step)+".dat","a")
                    file.write((2*"%.5E\t"+"\n")%(tracers[j][0], tracers[j][1]))    
                    file.close()

                    # 0) x-position                     \
                    # 1) y-position                      \
                    # 2) cell where the tracer is         --- Essential for advection of tracers ---
                    # 3) rank before advection           /
                    # 4) rank after advection           /

                    # --- Quantities to be advected ---
                    # 5) tau_xx
                    # 6) tau_xz
                    # 7) plastic strain
                    # 8) ocean material 
                    # 9) surface material
                    # 10) original depth of the tracer
                    # 11) composition (salty layer)
                    # 12) melt fraction
                    # 13) tracer original (0) or added (1)

            except:
                moving_tracers.append(tracers[j])
                tracers_in_cells[tracers[j][2]].remove(j)
                vacancy.append(j)
                tracers[j] = []

def find_tracers(moving_tracers, vacancy, tracers_to_find, name, method, tracers, t, step, dt, v, output_now):
    for k in range(0,size):
        for j in range (0,len(moving_tracers[k])):
            px = moving_tracers[k][j][0]
            py = moving_tracers[k][j][1]

            count1 = 0
            count2 = 0
            count3 = 0
            count4 = 0
            count5 = 0


            if (method=="Euler"):
                try:
                    v1 = v(Point(px,py))
                except:
                    v1 = [0.0,0.0]
                    count1+=1

                v1_x = MPI.sum(comm,v1[0])
                v1_y = MPI.sum(comm,v1[1])
                count1 = MPI.sum(comm,count1)     

                moving_tracers[k][j][0] = px + v1_y*float(dt)
                moving_tracers[k][j][1] = py + v1_y*float(dt)

            if (method=="RK2"):
                try:
                    v1 = v(Point(px,py))
                except:
                    v1 = [0.0,0.0]
                    count1+=1

                v1_x = MPI.sum(comm,v1[0])
                v1_y = MPI.sum(comm,v1[1])
                count1 = MPI.sum(comm,count1)

                try:
                    v2 = v(Point(px+v1_x*float(dt),py+v1_y*float(dt)))
                except:
                    v2 = [0.0,0.0]
                    count2+=1
                
                v2_x = MPI.sum(comm,v2[0])
                v2_y = MPI.sum(comm,v2[1])   
                count2 = MPI.sum(comm,count2)        

                moving_tracers[k][j][0] = px + (v1_x + v2_x)*float(dt)/2.0
                moving_tracers[k][j][1] = py + (v1_y + v2_y)*float(dt)/2.0

            if (method=="RK4"):
                try:
                    v1 = v(Point(px,py))
                except:
                    v1 = [0.0,0.0]
                    count1+=1

                v1_x = MPI.sum(comm,v1[0])
                v1_y = MPI.sum(comm,v1[1])
                count1 = MPI.sum(comm,count1)

                try:
                    v2 = v(Point(px+v1_x*float(dt)/2.0,py+v1_y*float(dt)/2.0))
                except:
                    v2 = [0.0,0.0]
                    count2+=1
                
                v2_x = MPI.sum(comm,v2[0])
                v2_y = MPI.sum(comm,v2[1])   
                count2 = MPI.sum(comm,count2)  

                try:
                    v3 = v(Point(px+v2_x*float(dt)/2.0,py+v2_y*float(dt)/2.0))
                except:
                    v3 = [0.0,0.0]
                    count3+=1

                v3_x = MPI.sum(comm,v3[0])
                v3_y = MPI.sum(comm,v3[1])
                count3 = MPI.sum(comm,count3)

                try:
                    v4 = v(Point(px+v3_x*float(dt),py+v3_y*float(dt)))
                except:
                    v4 = [0.0,0.0]
                    count4+=1

                v4_x = MPI.sum(comm,v4[0])             
                v4_y = MPI.sum(comm,v4[1])  
                count4 = MPI.sum(comm,count4)       

                moving_tracers[k][j][0] = px + (v1_x+2.0*v2_x+2.0*v3_x+v4_x)*float(dt)/6.0
                moving_tracers[k][j][1] = py + (v1_y+2.0*v2_y+2.0*v3_y+v4_y)*float(dt)/6.0

            # Find out the rank at the final (!!) position of the tracer, not earlier!
            try:
                v_test = v(Point(moving_tracers[k][j][0],moving_tracers[k][j][1]))
                moving_tracers[k][j][4] = rank
            except:
                moving_tracers[k][j][4] = 0
                count5+=1
            
            moving_tracers[k][j][4] = MPI.sum(comm,moving_tracers[k][j][4])
            count5 = MPI.sum(comm,count5) 

            not_found = False
            if (count1==size or count2==size or count3==size or count4==size or count5==size):
                # if count = size, the exeption has been performed by all processes, meaning that the tracer
                # left the domain, whatever the reason --> captures escape by both top and side boundaries
                not_found = True

            if (rank==int(moving_tracers[k][j][4]) and not_found==False):
                if (len(vacancy)>0): # if possible, put them in the vacancies
                    tracers[vacancy[0]] = moving_tracers[k][j]
                    tracers_to_find.append(vacancy[0]) 
                    del vacancy[0]

                else: # if no vacancy is available, put them at the end of the list
                    tracers.append(moving_tracers[k][j])
                    tracers_to_find.append(len(tracers)-1) # posledni pridany tracer

                # Output surface material tracers in every step
                # if (moving_tracers[k][j][9] > 0.0):
                #     file = open("data_"+name+"/tracers/surface_material/step_"+str(step)+".dat","a")
                #     file.write((2*"%.5E\t"+"%d"+"\n")%(moving_tracers[k][j][0],moving_tracers[k][j][1],moving_tracers[k][j][13]))    
                #     file.close()

                # if (moving_tracers[k][j][8] > 0.0):
                #     file = open("data_"+name+"/tracers/ocean_material/step_"+str(step)+".dat","a")
                #     file.write((2*"%.5E\t"+"\n")%(moving_tracers[k][j][0],moving_tracers[k][j][1]))    
                #     file.close()

                # Output all tracers on dedicated steps
                if (save_tracers==True and output_now==True):
                    file = open("data_"+name+"/tracers/step_"+str(step)+".dat","a")
                    # if (moving_tracers[k][j][11][0] == 1.0):
                    file.write((2*"%.5E\t"+"\n")%(moving_tracers[k][j][0],\
                                                moving_tracers[k][j][1]))    
                    file.close()



    
from dolfin import *
from m_parameters import *
from m_constants import *
from m_postproc import *
from pathlib import Path
import os 

class SaveFiles:
    def __init__(self, MeshClass, ElemClass):

        self.name = name
        if (Path("data_" + name).is_dir() == True):
            if (protect_directory == False):
                    if (rank == 0):
                        os.system("rm data_" + name + "/HDF5/*")
                        os.system("rm data_" + name + "/paraview/*.xdmf")
                        os.system("rm data_" + name + "/paraview/*.h5")
                        os.system("rm data_" + name + "/source_code/*.py")
                        os.system("rm data_" + name + "/boundary_data/*.dat")
                        os.system("rm data_" + name + "/tracers/*.dat")
            else:
                self.name += "_new"
                if (rank == 0):
                    print("\nWarning:\nName collision detected. New directory name is", str("data_" + self.name), ".\n")

        MPI.barrier(comm)
        directories = ["HDF5", "tracers", "paraview", "source_code"]
        for dir in directories:
            Path("data_" + self.name + "/" + dir).mkdir(parents = True, exist_ok = True)

        # --- Save the source code ---
        if (rank == 0):
            os.system("cp  main.py data_" + name + "/source_code")
            os.system("cp  m_*.py data_" + name + "/source_code")
            
        self.sDG0 = ElemClass.sDG0
        self.Temp = ElemClass.Temp
        self.v_k = ElemClass.v_k
        self.p_k = ElemClass.p_k
        self.comp_0 = ElemClass.comp_0
        self.comp_1 = ElemClass.comp_1
        self.comp_2 = ElemClass.comp_2
        self.heating = ElemClass.heating
        self.visc = ElemClass.visc
        self.log10_visc = ElemClass.log10_visc
        self.number_of_tracers = ElemClass.number_of_tracers
        self.mesh_ranks = ElemClass.mesh_ranks
        self.xm = ElemClass.xm
        self.u_mesh = ElemClass.u_mesh
        self.v_mesh = ElemClass.v_mesh
        self.h_bot = ElemClass.h_bot
        self.h_top = ElemClass.h_top
        self.n_bot = ElemClass.n_bot
        self.strain_rate_inv = ElemClass.strain_rate_inv
        self.iteration_error = ElemClass.iteration_error
        self.plastic_strain = ElemClass.plastic_strain
        self.yield_stress = ElemClass.yield_stress
        self.stress_dev_inv = ElemClass.stress_dev_inv
        self.stress_dev_inv_k = ElemClass.stress_dev_inv_k
        self.yield_function = ElemClass.yield_function
        self.eta_v = ElemClass.eta_v
        self.density = ElemClass.density
        self.z_function = ElemClass.z_function
        self.shear_modulus = ElemClass.shear_modulus

        self.mesh = MeshClass.mesh
        self.boundary_parts = MeshClass.boundary_parts
        self.save_mesh()

        self.Paraview_Dict = {}
        self.Function_Dict = {}

        # --- Initialize the text files ---
        file = open("data_" + self.name + "/statistics.dat","w")

        if (time_units == 1.0):
            file.write((2*"%s")%("Time (-)\t\t", "Step\t"))
        else:
            file.write((2*"%s")%("Time (" + time_units_string + ")\t\t", "Step\t"))

        for arg in stat_header:
            file.write(("%s\t\t")%(arg))
            
        file.write("\n")
        file.close()

        file = open("data_" + self.name + "/empty_cells.dat", "w")
        file.close()

        # --- Initialize the HDF5 files ---
        self.data_file = {}
        self.data_file['data'] = HDF5File(comm, "data_" + self.name + "/HDF5/data.h5", "w")

        stat_hdf5 = open("data_" + self.name + "/HDF5/data_timestamp.dat", 'w')
        stat_hdf5.write((4*"%s\t"+"\n")%("# HDF5", "step", "time step\t", "time"))    
        stat_hdf5.close()

        # --- Initialize Paraview files for initial condition and fill dictionaries ---
        for i in range(len(Paraview_Output_Ini)):
            self.Paraview_Dict[Paraview_Output_Ini[i]] = XDMFFile(comm,"data_" + self.name + "/paraview/initial_condition/" + Paraview_Output_Ini[i] + ".xdmf")
            self.Paraview_Dict[Paraview_Output_Ini[i]].parameters["flush_output"] = True
            self.Paraview_Dict[Paraview_Output_Ini[i]].parameters["rewrite_function_mesh"] = True

            if (Paraview_Output_Ini[i] == "temperature"):
                self.Function_Dict[Paraview_Output_Ini[i]] = self.Temp

        # --- Initialize Paraview files and fill dictionaries ---
        for i in range(len(Paraview_Output)):
            self.Paraview_Dict[Paraview_Output[i]] = XDMFFile(comm,"data_" + self.name + "/paraview/" + Paraview_Output[i] + ".xdmf")
            self.Paraview_Dict[Paraview_Output[i]].parameters["flush_output"] = True
            self.Paraview_Dict[Paraview_Output[i]].parameters["rewrite_function_mesh"] = True

            if (Paraview_Output[i] == "temperature"):
                self.Function_Dict[Paraview_Output[i]] = self.Temp
            
            if (Paraview_Output[i] == "velocity"):
                self.Function_Dict[Paraview_Output[i]] = self.v_k

            if (Paraview_Output[i] == "pressure"):
                self.Function_Dict[Paraview_Output[i]] = self.p_k

            if (Paraview_Output[i] == "strain_rate_inv"):
                self.Function_Dict[Paraview_Output[i]] = self.strain_rate_inv

            if (Paraview_Output[i] == "composition_0"):
                self.Function_Dict[Paraview_Output[i]] = self.comp_0

            if (Paraview_Output[i] == "composition_1"):
                self.Function_Dict[Paraview_Output[i]] = self.comp_1

            if (Paraview_Output[i] == "composition_2"):
                self.Function_Dict[Paraview_Output[i]] = self.comp_2

            if (Paraview_Output[i] == "tracers"):
                self.Function_Dict[Paraview_Output[i]] = self.number_of_tracers

            if (Paraview_Output[i] == "ranks"):
                self.Function_Dict[Paraview_Output[i]] = self.mesh_ranks

            if (Paraview_Output[i] == "viscosity_log"):
                self.Function_Dict[Paraview_Output[i]] = self.log10_visc

            if (Paraview_Output[i] == "viscosity"):
                self.Function_Dict[Paraview_Output[i]] = self.visc
            
            if (Paraview_Output[i] == "tidal_heating"):
                self.Function_Dict[Paraview_Output[i]] = self.heating
            
            if (Paraview_Output[i] == "melt_fraction"):
                self.Function_Dict[Paraview_Output[i]] = self.xm

            if (Paraview_Output[i] == "mesh_displacement"):
                self.Function_Dict[Paraview_Output[i]] = self.u_mesh

            if (Paraview_Output[i] == "mesh_velocity"):
                self.Function_Dict[Paraview_Output[i]] = self.v_mesh

            if (Paraview_Output[i] == "topography_bottom"):
                self.Function_Dict[Paraview_Output[i]] = self.h_bot

            if (Paraview_Output[i] == "topography_top"):
                self.Function_Dict[Paraview_Output[i]] = self.h_top
            
            if (Paraview_Output[i] == "normal_bottom"):
                self.Function_Dict[Paraview_Output[i]] = self.n_bot
            
            if (Paraview_Output[i] == "iteration_error"):
                self.Function_Dict[Paraview_Output[i]] = self.iteration_error

            if (Paraview_Output[i] == "plastic_strain"):
                self.Function_Dict[Paraview_Output[i]] = self.plastic_strain

            if (Paraview_Output[i] == "yield_stress"):
                self.Function_Dict[Paraview_Output[i]] = self.yield_stress

            if (Paraview_Output[i] == "stress_dev_inv"):
                self.Function_Dict[Paraview_Output[i]] = self.stress_dev_inv
            
            if (Paraview_Output[i] == "stress_dev_inv_k"):
                self.Function_Dict[Paraview_Output[i]] = self.stress_dev_inv_k

            if (Paraview_Output[i] == "yield_function"):
                self.Function_Dict[Paraview_Output[i]] = self.yield_function

            if (Paraview_Output[i] == "eta_v"):
                self.Function_Dict[Paraview_Output[i]] = self.eta_v

            if (Paraview_Output[i] == "density"):
                self.Function_Dict[Paraview_Output[i]] = self.density

            if (Paraview_Output[i] == "z_function"):
                self.Function_Dict[Paraview_Output[i]] = self.z_function

            if (Paraview_Output[i] == "shear_modulus"):
                self.Function_Dict[Paraview_Output[i]] = self.shear_modulus
    
    def write_statistic(self, t, step, stat_output, **kwargs):
        
        file = open("data_" + self.name + "/statistics.dat", "a")
        if (rank == 0):
            file.write(("%.5E\t\t" + "%d\t\t")%(float(t)/time_units, step))

        for arg in stat_output:
            if (arg == "nusselt"):
                value = nusselt(kwargs["q_cond_top"], kwargs["q_top"])
                if (rank == 0):
                    file.write(("%.5E\t\t")%(value))

            if (arg == "q_top"):
                value = kwargs["q_top"]
                if (rank == 0):
                    file.write(("%.5E\t\t")%(value))

            if (arg == "vrms"):
                value = rms_vel(kwargs["v"])
                if (rank == 0):
                    file.write(("%.5E\t\t")%(value))

            if (arg == "avg_h_bot"):
                if (rank == 0):
                    file.write(("%.5E\t\t")%(kwargs[arg]))

            if (arg == "h_top_max"):
                value = MPI.max(self.mesh.mpi_comm(), kwargs[arg].vector().max()) 
                if (rank == 0):
                    file.write(("%.5E\t\t")%(value))

            if (arg == "time"):
                value = kwargs["time"]
                if (rank == 0):
                    file.write(("%.5E\t\t")%(value))
            
            if (arg == "timestep"):
                if (rank == 0):
                    file.write(("%.5E\t\t")%(kwargs["timestep"]))
                
        if (rank == 0):
            file.write("\n")
        file.close()
        
    def Save_Paraview(self, t):
        for i in range(len(Paraview_Output)):
            self.Function_Dict[Paraview_Output[i]].rename(Paraview_Output[i], "")
            self.Paraview_Dict[Paraview_Output[i]].write(self.Function_Dict[Paraview_Output[i]], float(t)/time_units)

    def Save_Paraview_Ini(self):
        for i in range(len(Paraview_Output_Ini)):
            self.Function_Dict[Paraview_Output_Ini[i]].rename(Paraview_Output_Ini[i], "")
            self.Paraview_Dict[Paraview_Output_Ini[i]].write(self.Function_Dict[Paraview_Output_Ini[i]])

    def Save_HDF5(self, step_output, step, dt, t):
        for i in range(len(Paraview_Output)):
            self.data_file['data'].write(self.Function_Dict[Paraview_Output[i]], "/" + Paraview_Output[i], step) 
        
        self.data_file['data'].flush()

        if (rank == 0):
            stat_hdf5 = open("data_" + self.name + "/HDF5/data_timestamp.dat", 'a')
            stat_hdf5.write((2*"%d\t\t"+2*"%.5E\t\t"+"\n")%(step_output - 1, step, float(dt)/time_units, float(t)/time_units))
            stat_hdf5.close()

        # --- If the mesh is moving, save also the mesh ---
        # if (BC_vel_top == "free_surface" or BC_vel_bot == "free_surface"):
        self.Save_Mesh(step_output)

    def Load_HDF5(self, t, dt, *args):
        files_hdf5_in = {}
        files_hdf5_in['data'] = HDF5File(comm, "data_" + reload_name + "/HDF5/data.h5", "r")

        i = 0
        for arg in args:
            dataset = "/" + reload_HDF5_functions[i] + "/vector_%d" % (reload_HDF5_step)
            files_hdf5_in['data'].read(arg, dataset)
            i += 1

        # --- Read the time and timestep ---
        file_time = open("data_" + reload_name + "/HDF5/data_timestamp.dat")
        lines = file_time.readlines()
        sline = (lines[reload_HDF5_step + 1]).split("\t\t")

        dt.assign(float(sline[2])*time_units)
        if (restart_time == True):
            t.assign(0.0)
        else:
            t.assign(float(sline[3])*time_units)

        file_time.close()

    def save_mesh(self):
        mesh_file = XDMFFile(comm, "data_" + self.name + "/HDF5/mesh.xdmf")
        mesh_file.write(self.mesh)
        mesh_file.close()

        file = File("data_" + self.name + "/HDF5/subdomains.pvd")
        file.write(self.boundary_parts)

    def Save_Mesh(self, step_output):
        mesh_file = HDF5File(comm, "data_" + self.name + "/HDF5/meshes/mesh_" + str(int(step_output - 1)) + ".h5", "w")
        mesh_file.write(self.mesh, "/mesh")
        mesh_file.close()

    def Load_Mesh(self):
        self.mesh = Mesh()
        self.mesh_file = HDF5File(comm, "data_" + self.name + "/HDF5/meshes/mesh_" + str(int(reload_HDF5)) + ".h5", "r")
        self.mesh_file.read(mesh, "/mesh", True)
        self.mesh_file.close()
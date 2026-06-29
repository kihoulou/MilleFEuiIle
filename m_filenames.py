# --- Python modules ---
from dolfin import *
from pathlib import Path
import os 

# --- MilleFEuiIle modules ---
from m_parameters import *
from m_constants import *
from m_postproc import *

class SaveFiles:
    def __init__(self, MeshClass, ElemClass, use_tracers):
        self.name = name
        directories = ["HDF5",
                       "anim",
                       "img",
                       "topography",
                       "tracers",
                       "tracers/trajectories",
                       "paraview",
                       "paraview/initial_condition"]
        # /source_code created already in m_parameters.py

        if (Path("data_" + name).is_dir() == True):
            if (protect_directory == False):
                    if (rank == 0):
                        for dir in directories:
                            if (Path("data_" + name + "/" + dir).is_dir() == True):
                                os.system("rm -r data_" + name + "/" + dir)
                            else:
                                pass
                        
            else:
                self.name += "_new"
                if (rank == 0):
                    print("\nWarning:\nName collision detected. New directory name is", str("data_" + self.name), ".\n")

        MPI.barrier(comm)
        
        for dir in directories:
            Path("data_" + self.name + "/" + dir).mkdir(parents = True, exist_ok = True)

        # --- Save the source code ---
        if (rank == 0):
            os.system("cp  main.py data_" + name + "/source_code")

            files = ["boundary_conditions",
                     "check",
                     "constants",
                     "elements",
                     "equations",
                     "filenames", 
                     "incompatibility",
                     "interpolation",
                     "material_properties",
                     "melting",
                     "mesh",
                     "postproc",
                     "rheology",
                     "timestep",
                     "tracers"]
            
            # m_parameters.py copied already in m_parameters.py
            
            for f in files:
                os.system("cp  m_" + f + ".py data_" + name + "/source_code")
            
        self.sDG0 = ElemClass.sDG0
        self.Temp = ElemClass.Temp
        self.v_k = ElemClass.v_k
        self.p_k = ElemClass.p_k
        self.comp_0 = ElemClass.comp_0
        self.comp_1 = ElemClass.comp_1
        self.comp_2 = ElemClass.comp_2
        self.heating = ElemClass.heating
        self.visc = ElemClass.visc
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
        self.heating_shear = ElemClass.heating_shear
        self.cohesion = ElemClass.cohesion
        self.mechanisms = ElemClass.mechanisms
        self.diff_coef = ElemClass.diff_coef

        self.mesh = MeshClass.mesh
        self.boundary_parts = MeshClass.boundary_parts
        self.save_mesh_XDMF()

        self.Paraview_Dict = {}
        self.Paraview_Dict_Ini = {}
        self.Function_Dict = {}

        # --- Initialize the text files ---
        file = open("data_" + self.name + "/statistics.dat","w")

        if (time_units_string == "-"):
            file.write((2*"%s")%("Time (-)\t\t", "Step\t"))
        else:
            file.write((2*"%s")%("Time (" + time_units_string + ")\t\t", "Step\t\t"))

        for arg in stat_header:
            file.write(("%s\t\t")%(arg))
            
        file.write("\n")
        file.close()

        if (use_tracers == True):
            file = open("data_" + self.name + "/empty_cells.dat", "w")
            file.close()

        # --- Initialize the HDF5 files ---
        self.data_file = {}
        self.data_file['data'] = HDF5File(comm, "data_" + self.name + "/HDF5/data.h5", "w")

        stat_hdf5 = open("data_" + self.name + "/HDF5/data_timestamp.dat", 'w')
        stat_hdf5.write((4*"%s\t"+"\n")%("# HDF5 step", "Step\t", "Time step (" + time_units_string + ")\t", "Time (" + time_units_string + ")"))    
        stat_hdf5.close()

        # --- Initialize Paraview files for initial condition and fill dictionaries ---
        # --- Initialize Paraview files and fill dictionaries ---
        for i in range(len(paraview_output_ini)):
            self.Paraview_Dict_Ini[paraview_output_ini[i]] = XDMFFile(comm,"data_" + self.name + "/paraview/initial_condition/" + paraview_output_ini[i] + ".xdmf")
            self.Paraview_Dict_Ini[paraview_output_ini[i]].parameters["flush_output"] = True
            self.Paraview_Dict_Ini[paraview_output_ini[i]].parameters["rewrite_function_mesh"] = True

            self.assign_dictionaries(paraview_output_ini, self.Function_Dict, i)

        # --- Initialize Paraview files and fill dictionaries ---
        for i in range(len(paraview_output)):
            self.Paraview_Dict[paraview_output[i]] = XDMFFile(comm,"data_" + self.name + "/paraview/" + paraview_output[i] + ".xdmf")
            self.Paraview_Dict[paraview_output[i]].parameters["flush_output"] = True
            self.Paraview_Dict[paraview_output[i]].parameters["rewrite_function_mesh"] = True

            self.assign_dictionaries(paraview_output, self.Function_Dict, i)

    def assign_dictionaries(self, paraview_output, Function_Dict, i):
        if (paraview_output[i] == "temperature"):
            Function_Dict[paraview_output[i]] = self.Temp
        
        if (paraview_output[i] == "velocity"):
            Function_Dict[paraview_output[i]] = self.v_k

        if (paraview_output[i] == "pressure"):
            Function_Dict[paraview_output[i]] = self.p_k

        if (paraview_output[i] == "strain_rate_inv"):
            Function_Dict[paraview_output[i]] = self.strain_rate_inv

        if (paraview_output[i] == "composition_0"):
            Function_Dict[paraview_output[i]] = self.comp_0

        if (paraview_output[i] == "composition_1"):
            Function_Dict[paraview_output[i]] = self.comp_1

        if (paraview_output[i] == "composition_2"):
            Function_Dict[paraview_output[i]] = self.comp_2

        if (paraview_output[i] == "tracers"):
            Function_Dict[paraview_output[i]] = self.number_of_tracers

        if (paraview_output[i] == "ranks"):
            Function_Dict[paraview_output[i]] = self.mesh_ranks

        if (paraview_output[i] == "viscosity"):
            Function_Dict[paraview_output[i]] = self.visc
        
        if (paraview_output[i] == "tidal_heating"):
            Function_Dict[paraview_output[i]] = self.heating
        
        if (paraview_output[i] == "melt_fraction"):
            Function_Dict[paraview_output[i]] = self.xm

        if (paraview_output[i] == "mesh_displacement"):
            Function_Dict[paraview_output[i]] = self.u_mesh

        if (paraview_output[i] == "mesh_velocity"):
            Function_Dict[paraview_output[i]] = self.v_mesh

        if (paraview_output[i] == "topography_bottom"):
            Function_Dict[paraview_output[i]] = self.h_bot

        if (paraview_output[i] == "topography_top"):
            Function_Dict[paraview_output[i]] = self.h_top
        
        if (paraview_output[i] == "normal_bottom"):
            Function_Dict[paraview_output[i]] = self.n_bot
        
        if (paraview_output[i] == "iteration_error"):
            Function_Dict[paraview_output[i]] = self.iteration_error

        if (paraview_output[i] == "plastic_strain"):
            Function_Dict[paraview_output[i]] = self.plastic_strain

        if (paraview_output[i] == "yield_stress"):
            Function_Dict[paraview_output[i]] = self.yield_stress

        if (paraview_output[i] == "stress_dev_inv"):
            Function_Dict[paraview_output[i]] = self.stress_dev_inv
        
        if (paraview_output[i] == "stress_dev_inv_k"):
            Function_Dict[paraview_output[i]] = self.stress_dev_inv_k

        if (paraview_output[i] == "yield_function"):
            Function_Dict[paraview_output[i]] = self.yield_function

        if (paraview_output[i] == "eta_v"):
            Function_Dict[paraview_output[i]] = self.eta_v

        if (paraview_output[i] == "density"):
            Function_Dict[paraview_output[i]] = self.density

        if (paraview_output[i] == "z_function"):
            Function_Dict[paraview_output[i]] = self.z_function

        if (paraview_output[i] == "shear_modulus"):
            Function_Dict[paraview_output[i]] = self.shear_modulus

        if (paraview_output[i] == "cohesion"):
            Function_Dict[paraview_output[i]] = self.cohesion

        if (paraview_output[i] == "mechanisms"):
            Function_Dict[paraview_output[i]] = self.mechanisms

        if (paraview_output[i] == "diff_coef"):
            Function_Dict[paraview_output[i]] = self.diff_coef
    

    def write_statistic(self, t, step, stat_output, time_units, **kwargs):
        
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
        
    def save_paraview(self, t, time_units):
        for i in range(len(paraview_output)):
            self.Function_Dict[paraview_output[i]].rename(paraview_output[i], "")
            self.Paraview_Dict[paraview_output[i]].write(self.Function_Dict[paraview_output[i]], float(t)/time_units)

    def save_paraview_ini(self):
        for i in range(len(paraview_output_ini)):
            self.Function_Dict[paraview_output_ini[i]].rename(paraview_output_ini[i], "")
            self.Paraview_Dict_Ini[paraview_output_ini[i]].write(self.Function_Dict[paraview_output_ini[i]])

    def save_HDF5(self, step_output, step, dt, t, time_units):
        for i in range(len(paraview_output)):
            self.data_file['data'].write(self.Function_Dict[paraview_output[i]], "/" + paraview_output[i], step) 
        
        self.data_file['data'].flush()

        if (rank == 0):
            stat_hdf5 = open("data_" + self.name + "/HDF5/data_timestamp.dat", 'a')
            stat_hdf5.write((2*"%d\t\t"+2*"%.5E\t\t"+"\n")%(step_output - 1, step, float(dt)/time_units, float(t)/time_units))
            stat_hdf5.close()

        # --- Always save the mesh as well - used for reloading ---
        self.save_mesh_HDF5(step_output)

    def load_HDF5(self, t, dt, time_units, name, *args):

        files_hdf5_in = {}
        files_hdf5_in['data'] = HDF5File(comm, reload_name + "/HDF5/data.h5", "r")

        file_time = open(reload_name + "/HDF5/data_timestamp.dat")
        lines = file_time.readlines()
        
        if (reload_step == "last"):
            step = len(lines) - 2
        else:
            step = reload_step

        i = 0
        for arg in args:

            dataset = "/" + name + "/vector_" + str(step)
            files_hdf5_in['data'].read(arg, dataset)
            i += 1

        # --- Read the time and timestep ---
        sline = (lines[step + 1]).split("\t\t") #  +1 because of the header
        dt.assign(float(sline[2])*time_units)
        if (restart_time == True):
            t.assign(0.0)
        else:
            t.assign(float(sline[3])*time_units)

        file_time.close()

    def save_mesh_XDMF(self):
        mesh_file = XDMFFile(comm, "data_" + self.name + "/HDF5/mesh.xdmf")
        mesh_file.write(self.mesh)
        mesh_file.close()

        file = File("data_" + self.name + "/HDF5/subdomains.pvd")
        file.write(self.boundary_parts)

    def save_mesh_HDF5(self, step_output):
        mesh_file = HDF5File(comm, "data_" + self.name + "/HDF5/meshes/mesh_" + str(int(step_output - 1)) + ".h5", "w")
        mesh_file.write(self.mesh, "/mesh")
        mesh_file.write(self.boundary_parts, "/subdomains")
        mesh_file.close()
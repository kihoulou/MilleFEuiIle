from m_parameters import *
from m_constants import *
from m_postproc import *
from pathlib import Path

class SaveFiles:
    def __init__(self, MeshClass, EqClass):

        directories = ["HDF5", "tracers", "paraview", "source_code"]
        for dir in directories:
            Path("data_" + name + "/" + dir).mkdir(parents = True, exist_ok = True)

        self.Temp = EqClass.Temp
        self.v_k = EqClass.v_k
        self.comp_0 = EqClass.comp_0
        self.comp_1 = EqClass.comp_1
        self.comp_2 = EqClass.comp_2
        self.number_of_tracers = EqClass.number_of_tracers
        self.mesh_ranks = EqClass.mesh_ranks

        self.mesh = MeshClass.mesh
        self.Save_Mesh = MeshClass.Save_Mesh

        self.Paraview_Dict = {}
        self.Function_Dict = {}

        # --- Initialize the text files ---
        file = open("data_" + name + "/statistics.dat","w")

        if (time_units == 1.0):
            file.write((2*"%s")%("Time (-)\t\t", "Step\t"))
        else:
            file.write((2*"%s")%("Time ("+str(time_units)+")\t\t", "Step\t"))

        for arg in stat_header:
            file.write(("%s\t\t")%(arg))
            
        file.write("\n")
        file.close()

        file = open("data_" + name + "/empty_cells.dat", "w")
        file.close()

        # --- Initialize the HDF5 files ---
        self.data_file = {}
        self.data_file['data'] = HDF5File(comm, "data_" + name + "/HDF5/data.h5", "w")

        stat_hdf5 = open("data_"+name+"/HDF5/data_timestamp.dat", 'w')
        stat_hdf5.write((4*"%s\t"+"\n")%("# HDF5", "step", "time step\t", "time"))    
        stat_hdf5.close()

        # --- Initialize Paraview files and fill dictionaries ---
        for i in range(len(Paraview_Output)):
            self.Paraview_Dict[Paraview_Output[i]] = XDMFFile(comm,"data_" + name + "/paraview/" + Paraview_Output[i] + ".xdmf")
            self.Paraview_Dict[Paraview_Output[i]].parameters["flush_output"] = True
            self.Paraview_Dict[Paraview_Output[i]].parameters["rewrite_function_mesh"] = True

            if (Paraview_Output[i] == "temperature"):
                self.Function_Dict[Paraview_Output[i]] = self.Temp
            
            if (Paraview_Output[i] == "velocity"):
                self.Function_Dict[Paraview_Output[i]] = self.v_k

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
        
    
    def write_statistic(self, t, step, stat_output, **kwargs):
        
        file = open("data_" + name + "/statistics.dat", "a")
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

    def Save_HDF5(self, step_output, step, dt, t):
        for i in range(len(Paraview_Output)):
            self.data_file['data'].write(self.Function_Dict[Paraview_Output[i]], "/" + Paraview_Output[i], step) 
        
        self.data_file['data'].flush()

        if (rank == 0):
            stat_hdf5 = open("data_"+name+"/HDF5/data_timestamp.dat", 'a')
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
            dataset = "/" + Functions_Input[i] + "/vector_%d" % (reload_HDF5)
            files_hdf5_in['data'].read(arg, dataset)
            i += 1

        # --- Read the time and timestep ---
        file_time = open("data_" + reload_name + "/HDF5/data_timestamp.dat")
        lines = file_time.readlines()
        sline = (lines[reload_HDF5 + 1]).split("\t\t")

        dt.assign(float(sline[2])*time_units)
        if (restart_time == True):
            t.assign(0.0)
        else:
            t.assign(float(sline[3])*time_units)

        file_time.close()
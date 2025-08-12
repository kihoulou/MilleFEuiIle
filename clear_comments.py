import subprocess

subprocess.run(["rm  clean/*"], shell=True) 

files = ["main",
         "m_boundary_conditions",
         "m_check",
         "m_constants",
         "m_equations",
         "m_elements",
         "m_filenames",
         "m_interpolation",
         "m_incompatibility",
         "m_material_properties",
         "m_melting",
         "m_mesh",
         "m_parameters_europa",
         "m_postproc",
         "m_rheology",
         "m_timestep",
         "m_tracers"]

for file in files:

    file_in = open(file + ".py", "r")
    file_out = open("clean/" + file +".py", "w")

    comment = False

    while line := file_in.readline():

        if ('"""' in line):
            comment = not comment
            continue

        if (comment == False):
            file_out.write(line)

        else:
            pass

    print("File " + file + ".py finished.")
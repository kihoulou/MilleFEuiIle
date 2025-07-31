# Number of cores
n_cores=8

# Parameters file
param_file="m_parameters_tectonics_VEP.py"
cp $param_file m_parameters.py

mpirun -n $n_cores python main.py

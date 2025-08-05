# Number of cores
n_cores=4

# Parameters file
param_file="m_parameters_convection.py"
cp $param_file m_parameters.py

mpirun -n $n_cores python main.py

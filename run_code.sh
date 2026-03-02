# --- Number of cores ---
n_cores=8

# --- Unique extension to parameter and main file ---
ext="ganymede"

main_file1="main.py"
main_file2="main_"$ext".py"

param_file1="m_parameters.py"
param_file2="m_parameters_"$ext".py"

out_file1=$ext".out"
out_file2=$ext"_e.out"

cp $param_file2 $param_file1
cp $main_file1 $main_file2

# mpirun -n $n_cores python $main_file2 > $out_file1 2> $out_file2
mpirun -n $n_cores python $main_file2


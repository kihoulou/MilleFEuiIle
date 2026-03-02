# --- Number of cores ---
n_cores=1

# --- Unique extension to parameter and main file (m_parameters_*.py) ---
name="docs"

main_file1="main.py"
main_file2="main_"$name".py"

param_file1="m_parameters.py"
param_file2="m_parameters_"$name".py"

out_file1=$name".out"
out_file2=$name"_e.out"

cp $param_file2 $param_file1
cp $main_file1 $main_file2

if [ $n_cores -eq 1 ]; then
    python $main_file2 > $out_file1 2> $out_file2&
else
    mpirun -n $n_cores python $main_file2 > $out_file1 2> $out_file2&
fi



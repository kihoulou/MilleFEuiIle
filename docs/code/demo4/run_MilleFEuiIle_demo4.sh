#-------------------------------------------------
# Launching script for MilleFEuiIle version 1.0.1
#-------------------------------------------------

# --- Loop over parameters? ---
# loop=0 -> NO, the script launches the code once
# loop=1 -> YES, the script launches a parametric sweep with one simulation at a time, running at n_cores cores
# loop=2 -> YES, the script launches as much as possible until n_cores_max, launches next as soon the previous ones finish
loop=0

# --- Number of cores ---
n_cores=4

# --- Maximum number of cores the script can take ---
# relevant only if loop=2
n_cores_max=8

# --- Running on background or foreground ---
# background=0 -> running on foreground and output to screen
# background=1 -> running on background and output to file
# choice loop=1 or loop=2 below always write to file
background=0

# --- Unique extension to parameter and main file (m_parameters_*.py) ---
name="demo4_reload"

main_file1="main.py"
param_file1="m_parameters.py"

# --- Function run_MilleFEuiIle 
run_MilleFEuiIle() {
    main_file2="main_MF_"$name".py"
    param_file2="m_parameters_"$name".py"

    out_file1=$name".out"
    out_file2=$name"_e.out"

    cp $param_file2 $param_file1
    cp $main_file1 $main_file2

    if [ $n_cores -eq 1 ]; then
        if [ $background -eq 0 && $loop -eq 0 ]; then
            python $main_file2
        elif [ $background -eq 1 && $loop -eq 0 ]; then
            python $main_file2 > $out_file1 2> $out_file2&
        fi
    else
        if [ $background -eq 0 ] && [ $loop -eq 0 ]; then
            mpirun -n $n_cores python $main_file2
        elif [ $background -eq 1 ] && [ $loop -eq 0 ]; then
            mpirun -n $n_cores python $main_file2 > $out_file1 2> $out_file2&
        fi
    fi
}

run_MilleFEuiIle_loop() {
    main_file2="main_MF_"$name"_par"$par".py"
    param_file2="m_parameters_"$name".py"

    out_file1=$name"_par"$par".out"
    out_file2=$name"_par"$par"_e.out"

    cp $param_file2 $param_file1
    cp $main_file1 $main_file2

    if [ $loop -eq 1 ]; then
        if [ $n_cores -eq 1 ]; then
            python $main_file2 $par > $out_file1 2> $out_file2
        else
            mpirun -n $n_cores python $main_file2 $par > $out_file1 2> $out_file2
        fi

    elif [ $loop -eq 2 ]; then
        if [ $n_cores -eq 1 ]; then
            python $main_file2 $par > $out_file1 2> $out_file2&
        else
            mpirun -n $n_cores python $main_file2 $par > $out_file1 2> $out_file2&
        fi
    fi
}

# --- Running the code a single time---
if [ $loop -eq 0 ]; then
    run_MilleFEuiIle

# --- Running the code a loop over parameters (always at n_cores cores) ---
elif [ $loop -eq 1 ]; then
    for par in 1
        do
            run_MilleFEuiIle_loop par
        done

# --- Running the code a loop over parameters (up to n_cores_max cores) ---
elif [ $loop -eq 2 ]; then
    # --- Parameter loop ---
    for par in 13 15 14 16
        do
            # --- Waiting loop ---
            while true
                do 
                    # --- Check every second if there are enough free processors ---
                    sleep 1       
                    if [ $(ps -ef | grep -v grep | grep main_MF | wc -l) -lt $n_cores_max ]; then

                        # --- Launch once and exit the waiting loop ---
                        run_MilleFEuiIle_loop par
                        break
                    
                    fi
                done
        done
fi
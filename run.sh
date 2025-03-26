# Number of cores
n_cores=8

# Vertical resolution
z_div=100

# Tectonic mode
mode="ext"

logfile="simulations_"$mode"_res"$z_div".log"

echo -e "Loop over maximum grain size \n" > "$logfile"

echo -e "PID of this script:" $$ "\n" >> "$logfile"

start0=$(date +%s)

for i in {0..40..10}
    do        
        main_name="main_"$mode"_res"$z_div"_phi"$i".py"
        cp main.py $main_name
        echo "Running mpirun -n $n_cores $main_name $mode $z_div $i" >> "$logfile"
        
        start=$(date +%s)
        mpirun -n $n_cores python $main_name $mode $z_div $i > "$mode"_res"$z_div"_phi"$i".out 2> "$mode"_res"$z_div"_phi"$i"_e.out
        end=$(date +%s)
        
        rm $main_name
        
        ((hours=("$end"-"$start")/3600))
        ((min=(("$end"-"$start") - "$hours"*3600)/60))
        ((sec=(("$end"-"$start") - "$hours"*3600 - "$min"*60)))

        ((hours0=("$end"-"$start0")/3600))
        ((min0=(("$end"-"$start0") - "$hours0"*3600)/60))
        ((sec0=(("$end"-"$start0") - "$hours0"*3600 - "$min0"*60)))

        echo -e "Done! Elapsed time $hours h $min min $sec sec | Total time $hours0 h $min0 min $sec0 sec.\n" >> "$logfile"
done

echo -e "Script finished." >> "$logfile"

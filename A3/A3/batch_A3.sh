#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 60:00
#SBATCH -J quicksort

################################################################################
# Batch script for Assignment 3 in Parallel and Distributed Programming spring 2023.
# Note that the following commands need to be run before this script:
#   module load gcc openmpi
#   make clean
#   make

# Runs the program with different input on different number of processes to check
# strong scalability. (Weak scalability can be checked from the output of the
# script with some manual work.)
# Author: Jacob Malmenstedt 2023
################################################################################

num_proc=(
    1 
    2 
    4 
    8 
    16
)
pivot_method=(  
    1 
    2 
    3
)
data_dir=/proj/uppmax2023-2-13/nobackup/qsort_indata
output_file="output_to_remove.txt"
commands=(
	"./quicksort ${data_dir}/input125000000.txt ${output_file}"
	# "./quicksort ${data_dir}/input250000000.txt ${output_file}"
	# "./quicksort ${data_dir}/input500000000.txt ${output_file}"
	# "./quicksort ${data_dir}/input1000000000.txt  ${output_file}"
    # "./quicksort ${data_dir}/input2000000000.txt  ${output_file}"
)

for pivot in ${pivot_method[@]}; do
    echo "Checking strong scalability with pivot method $pivot"
    for p in ${num_proc[@]}; do
        for (( i = 0; i < ${#commands[@]}; i++ )); do
            echo "Running ${commands[$i]} with $p processes"
            mpirun --bind-to none -np $p ${commands[$i]} $pivot
        done
    done
done
echo "OK"




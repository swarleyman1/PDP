#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 1:00:00
#SBATCH -J quicksort_v2

################################################################################
# Batch script for Assignment 3 in Parallel and Distributed Programming spring 2023.
# Note that the following commands need to be run before this script:
#   module load gcc openmpi
#   make clean
#   make quicksort_v2

# Runs the program with different input on different number of processes to check
# strong scalability. (Weak scalability can be checked from the output of the
# script with some manual work.)
# Author: Jacob Malmenstedt 2023
################################################################################

num_proc=(1 2 4 8 16)

data_dir=/proj/uppmax2023-2-13/nobackup/qsort_indata
output_file="output_to_remove.txt"
commands=(
	# "./quicksort_v2 ${data_dir}/input125000000.txt ${output_file} 1"
	"./quicksort_v2 ${data_dir}/input250000000.txt ${output_file} 1"
	# "./quicksort_v2 ${data_dir}/input500000000.txt ${output_file} 1"
	# "./quicksort_v2 ${data_dir}/input1000000000.txt  ${output_file} 1"
    # "./quicksort_v2 ${data_dir}/input2000000000.txt  ${output_file} 1"
)

for p in ${num_proc[@]}; do
    for (( i = 0; i < ${#commands[@]}; i++ )); do
        echo "Running ${commands[$i]} with $p processes"
        mpirun --bind-to none -np $p ${commands[$i]}
    done
done
echo "OK"




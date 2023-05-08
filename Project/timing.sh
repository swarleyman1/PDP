#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 20:00
#SBATCH -J mc

################################################################################
# Batch script for Project in Parallel and Distributed Programming spring 2023.

# For running the code:
# mpirun -np <number of processes> ./mc <number of simulations per process> <output file>

# Runs the program with different number of simulations and processes to check
# scalability.
################################################################################

num_proc=(1 2 4 8 16)
num_sim=(1000000 2000000 4000000 8000000 16000000)
output_file="output_to_remove.txt"

for proc in ${num_proc[@]}; do
    for sim in ${num_sim[@]}; do
        echo "Running with $proc processes and $sim simulations"
        mpirun --bind-to none -np $proc ./mc $sim $output_file
    done
done
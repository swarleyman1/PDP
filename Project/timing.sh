#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 3:00:00
#SBATCH -J mc

################################################################################
# Batch script for Project in Parallel and Distributed Programming spring 2023.

# For running the code:
# mpirun -np <number of processes> ./mc <number of simulations per process> <output file>

# Runs the program with different number of simulations and processes to check
# scalability.
################################################################################


num_proc=(16 8 4 2)
num_sim=(100000 200000 400000 800000)
output_file="output_to_remove.txt"

for proc in "${num_proc[@]}"; do
    for sim in "${num_sim[@]}"; do
        sim_per_proc=$((sim / proc))  # Calculate number of simulations per process
        echo "Running with $proc processes and $sim_per_proc simulations per process"
        echo "In total $sim simulations"
        mpirun --bind-to none -np "$proc" ./mc "$sim_per_proc" "$output_file"
    done
done

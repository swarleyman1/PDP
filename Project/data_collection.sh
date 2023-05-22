#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 2:00:00
#SBATCH -J mc

echo "Running with 16 processes and 1000000 simulations"
mpirun --bind-to none -np 16 ./mc 62500 data1000000.csv
echo "Running with 16 processes and 2000000 simulations"
mpirun --bind-to none -np 16 ./mc 125000 data2000000.csv
echo "Running with 16 processes and 4000000 simulations"
mpirun --bind-to none -np 16 ./mc 250000 data4000000.csv
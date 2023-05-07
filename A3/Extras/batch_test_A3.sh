#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 10:00
#SBATCH -J quicksort

echo "Running quicksort_v2"
echo "running with pivot method 1"
mpirun --bind-to none -np 8 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output_to_remove.txt 1

echo "running with pivot method 2"
mpirun --bind-to none -np 8 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output_to_remove.txt 2

echo "running with pivot method 3"
mpirun --bind-to none -np 8 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output_to_remove.txt 3

echo "All done!"



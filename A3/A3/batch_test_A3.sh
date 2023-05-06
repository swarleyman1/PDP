#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 5:00
#SBATCH -J quicksort


mpirun -np 2 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input2000000000.txt output_to_remove.txt 3

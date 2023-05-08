# Monte Carlo SSA (Stochastic Simulation Algorithm)

This is an implementation of the Monte Carlo SSA with OpenMPI multiprocessing for the course Parallel and Distributed Programming at Uppsala University. 

The code will write the dats for a histogram of the final number of susceptible humans to the specified output file after all experiments are run.

Written by Jacob Malmenstedt in 2023. 

## Usage

The code takes two input parameters:

- `n`: the number of simulations to run in each process
- `output file`: The file to which the results should be written to

## Propensity Function

The propensity function used in this implementation is set up to simulate a malaria outbreak. However, this can be changed in the `prop.c` code file to simulate other stochastic processes.

## Running

To run the code, run the following command:

```bash
mpirun -np <number of processes> ./mc <number of simulations> <output file>
```


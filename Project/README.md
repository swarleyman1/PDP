# Monte Carlo SSA (Stochastic Simulation Algorithm)

This is an implementation of the Monte Carlo SSA with OpenMPI multiprocessing for the course Parallel and Distributed Programming at Uppsala University

Written by Jacob Malmenstedt in 2023. 

## Usage

The code takes one input parameter:

- `n`: the number of simulations to run in each process


## Propensity Function

The propensity function used in this implementation is set up to simulate a malaria outbreak. However, this can be changed in the `prop.c` code file to simulate other stochastic processes.

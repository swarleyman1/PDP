#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "prop.h"
#include "update_state.h"

/****************************************************************************
 * Code for the project in the course PDP at Uppsala University 2023        *
 *                                                                          *
 * An implementation of the Monte Carlo SSA with OpenMPI multiprocessing    *
 * Author: Jacob Malmenstedt                                                *
 *                                                                          *
 ****************************************************************************/

// Upgrades to be implemented:
// Use openMP to parallelize the execution of the experiments on each process.


/* For higher grades:
All processes work fully in parallel, however, as the time steps are chosen randomly, the processes
proceed to reach the final time asynchronously. 
Record the wall clock time after passing time 25,50,75,100 (per process)
so that at the end you can output the average time per processor for
each time sub-interval. You are encouraged to use the MPI one-sided put/get functionality.*/

// Constants
#define OUTPUT 1    // Flag for output to file or not
#define R 15        // Number of reactions
#define Q 7         // Number of quantities
#define T 100       // Maximum time

// Function declarations
void prop(int *x, double *w);
static inline double rand_zero_to_one();
static inline int select_reaction(double *q, double a0, int length);
static int gillespieSSA(int *x, double *w, double *q, const int *x0);

// Main function
int main(int argc, char *argv[])
{

    // Check if correct number of arguments
    #if OUTPUT
    if (argc != 3)
    {
        printf("Usage: %s <number of iterations> <output file>\n", argv[0]);
        return 1;
    }
    #else
    if (argc != 2)
    {
        printf("Usage: %s <number of iterations>\n", argv[0]);
        return 1;
    }
    #endif

    // Read command line arguments
    int N = atoi(argv[1]);    // Number of iterations
    #if OUTPUT
    char *filename = argv[2]; // Output file name
    #endif
    
    // Initialize MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    //MPI_Status status;
    //MPI_Request send_request, recv_request;
    MPI_Comm_size(MPI_COMM_WORLD, &size);   // Number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // Rank of current process

    // Initialize
    const int n = N / size;                             // Number of iterations per process
    const int x0[] = {900, 900, 30, 330, 50, 270, 20};  // Initial state vector
    int *local_X = (int *)malloc(n * sizeof(int));      // Local output vector
    int *x = (int *)malloc(Q * sizeof(int));            // State vector
    double *w = (double *)malloc(R * sizeof(double));   // Propensity vector
    double *q = (double *)malloc(R * sizeof(double));   // Cumulative propensity vector

    // Initialize random number generator
    srand(time(NULL) + rank);

    // Start timer
    double time = MPI_Wtime();



    // Run experiments n times on each process
    for (int i = 0; i < n; i++)
    {
        local_X[i] = gillespieSSA(x, w, q, x0);
    }



    // Find max and min element in local output vector
    int local_max = INT_MIN;
    int local_min = INT_MAX;
    for (int i = 0; i < n; i++)
    {
        if (local_X[i] > local_max)
        {
            local_max = local_X[i];
        }
        if (local_X[i] < local_min)
        {
            local_min = local_X[i];
        }
    }

    // Find global max and min element
    int global_max, global_min;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Create 20 histogram bins
    int *bins = (int *)malloc(20 * sizeof(int));
    double bin_size = (global_max - global_min) / 20.0;
    for (int i = 0; i < 20; i++)
    {
        bins[i] = global_min + i * bin_size;
    }

    // Count number of elements in each bin
    int *local_bin_counts = (int *)calloc(20, sizeof(int));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 20; j++)
        {
            if (local_X[i] >= bins[j] && local_X[i] < bins[j + 1])
            {
                local_bin_counts[j]++;
            }
        }
        if (local_X[i] == global_max)
        {
            local_bin_counts[19]++;
        }
    }

    // Find global bin counts
    int *global_bin_counts = (int *)calloc(20, sizeof(int));
    MPI_Reduce(local_bin_counts, global_bin_counts, 20, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Stop timer
    time = MPI_Wtime() - time;

    // Print max time
    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("%f\n", max_time);
    }

    // Write to file
    #if OUTPUT
    if (rank == 0)
    {
        FILE *fp;
        fp = fopen(filename, "w");
        fprintf(fp, "Nr iters: \t %d\n", N);
        fprintf(fp, "Max time: \t %.4f \n", max_time);
        fprintf(fp, "Min: \t\t %d\n", global_min);
        fprintf(fp, "Max: \t\t %d\n \n", global_max);
        fprintf(fp, "Bin \t\t Count\n");
        fprintf(fp, "------------------\n");
        for (int i = 0; i < 19; i++)
        {
            fprintf(fp, "%d - %d \t %d\n", bins[i], bins[i + 1], global_bin_counts[i]);
        }
        fprintf(fp, "%d - %d \t %d\n", bins[19], global_max, global_bin_counts[19]);
        fclose(fp);
    }
    #endif

    // Free memory
    free(local_X);
    free(bins);
    free(local_bin_counts);
    free(global_bin_counts);
    free(x);
    free(w);
    free(q);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}

// Function definitions
static inline double rand_zero_to_one()
{
    // Generate random double between 0 and 1
    return (double)rand() / (double)RAND_MAX;
}

static inline int select_reaction(double *q, double a0, int length)
{
    // Multiply a0 with random double between 0 and 1
    double u2a0 = rand_zero_to_one() * a0;

    // Find index of first element in q that is larger than u2*a0
    int i = 0;
    while (i < length && q[i] < u2a0)
    {
        i++;
    }

    return i;
}

static int gillespieSSA(int *x, double *w, double *q, const int *x0)
{
    double t = 0.0;  // Time
    double tau;      // Time step
    double a0;       // Total propensity
    int reaction;    // Reaction index

    memcpy(x, x0, Q * sizeof(int));             // Copy initial state to state vector


    // Main loop
    while (t < T)
    {
        // Calculate propensities
        prop(x, w);

        // Calculate cumulative propensities
        q[0] = w[0];
        for (int i = 1; i < R; i++)
        {
            q[i] = q[i - 1] + w[i];
        }

        // Calculate total propensity
        a0 = q[R - 1];

        // Calculate time step
        tau = -log(rand_zero_to_one()) / a0;

        // Select reaction
        reaction = select_reaction(q, a0, R);

        // Update state
        update_state(x, reaction);

        // Update time
        t += tau;
    }

    return x[0];
}

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
 * (See README.md for more information)                                     *
 *                                                                          *
 * Author: Jacob Malmenstedt                                                *
 *                                                                          *
 ****************************************************************************/

// Constants
#define OUTPUT 2 // Flag for output: 0 = no output, 1 = human readable, 2 = csv
#define R 15     // Number of reactions
#define Q 7      // Number of quantities
#define T 100    // Maximum time

// Function declarations
void prop(int *x, double *w);
static inline double rand_zero_to_one();
static inline int select_reaction(double *q, double a0, int length);
static int gillespieSSA(int *x, double *w, double *q, const int *x0, const int *checkpoints, const int rank, const int size, double *timings);

// Main function
int main(int argc, char *argv[])
{

// Check if correct number of arguments
#if OUTPUT > 0
    if (argc != 3)
    {
        printf("Usage: %s <number of iterations> <output file>\n", argv[0]);
        return 1;
    }
#else
    if (argc < 2)
    {
        printf("Usage: %s <number of iterations>\n", argv[0]);
        return 1;
    }
#endif

    // Read command line arguments
    int n = atoi(argv[1]); // Number of iterations per process
#if OUTPUT > 0
    char *filename = argv[2]; // Output file name
#endif

    // Initialize MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Rank of current process

    // Initialize
    const int N = n * size;                                                      // Number of iterations in total
    const int checkpoints[4] = {25, 50, 75, 100};                                // Checkpoints when to record runtimes
    const size_t num_checkpoints = sizeof(checkpoints) / sizeof(checkpoints[0]); // Number of times to record
    double runtimes[num_checkpoints * size];                                     // Array to store runtimes
    double local_runtimes[num_checkpoints];                                      // Array to store local runtimes
    const int x0[] = {900, 900, 30, 330, 50, 270, 20};                           // Initial state vector
    int local_X[n];                                                              // Local output vector
    int x[Q];                                                                    // State vector
    double w[R];                                                                 // Propensity vector
    double q[R];                                                                 // Cumulative propensity vector

    // Initialize arrays with 0.0
    memset(runtimes, 0.0, sizeof(runtimes));
    memset(local_runtimes, 0.0, sizeof(local_runtimes));

    // Initialize random number generator
    srand(time(NULL) + rank);

    // Start timer
    double start_time = MPI_Wtime();

    // Run experiments n times on each process
    for (int i = 0; i < n; i++)
    {
        local_X[i] = gillespieSSA(x, w, q, x0, checkpoints, rank, size, local_runtimes);
    }

    double mid_time = MPI_Wtime();

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
    int bins[20];
    const double bin_size = (global_max - global_min) / 20.0;
    for (int i = 0; i < 20; i++)
    {
        bins[i] = global_min + i * bin_size;
    }

    // Count number of elements in each bin
    int local_bin_counts[20];
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
    int global_bin_counts[20];
    MPI_Reduce(local_bin_counts, global_bin_counts, 20, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Stop timer
    double end_time = MPI_Wtime();

    // Calculate runtimes
    const double tot_time = end_time - start_time;
    const double comm_time = end_time - mid_time;

    // Reduce runtimes
    double avg_time, max_time, min_time;
    double avg_comm_time, max_comm_time, min_comm_time;
    MPI_Reduce(&tot_time, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tot_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tot_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    MPI_Reduce(&comm_time, &avg_comm_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time, &min_comm_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    // Gather checkpoint timings
    MPI_Gather(local_runtimes, num_checkpoints, MPI_DOUBLE, runtimes, num_checkpoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Print time statistics
    if (rank == 0)
    {
        printf("total time for %d iterations:\n", N);
        printf("max: %.4f, min: %.4f, avg: %.4f\n\n", max_time, min_time, avg_time / size);
        printf("communication time for %d iterations:\n", N);
        printf("max: %.4f, min: %.4f, avg: %.4f\n\n", max_comm_time, min_comm_time, avg_comm_time / size);
        printf("RANK\t\t");
        for (int i = 0; i < size; i++)
        {
            printf("%d         ", i);
        }
        printf("\n");
        for (int i = 0; i < num_checkpoints; i++)
        {
            printf("Checkpoint %d\t", checkpoints[i]);
            for (int j = 0; j < size; j++)
            {
                printf("%.4fms  ", runtimes[j * num_checkpoints + i] * 1000 / N);
            }
            printf("\n");
        }
    }

// Write to file
#if OUTPUT == 1 // Human readable output
    if (rank == 0)
    {
        FILE *fp;
        fp = fopen(filename, "w");
        fprintf(fp, "Nr iters: \t %d\n", N);
        fprintf(fp, "Nr procs: \t %d\n", size);
        fprintf(fp, "Avg time: \t %.4f \n", avg_time / size);
        fprintf(fp, "Max time: \t %.4f \n", max_time);
        fprintf(fp, "Min time: \t %.4f \n", min_time);
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
#elif OUTPUT == 2 // CSV output
    if (rank == 0)
    {
        FILE *fp;
        fp = fopen(filename, "w");
        fprintf(fp, "%d, %d, %.4f, %.4f, %.4f, %d, %d\n", N, size, avg_time / size, max_time, min_time, global_min, global_max);
        for (int i = 0; i < 19; i++)
        {
            fprintf(fp, "%d, %d, %d\n", bins[i], bins[i + 1], global_bin_counts[i]);
        }
        fprintf(fp, "%d, %d, %d\n", bins[19], global_max, global_bin_counts[19]);
        fclose(fp);
    }
#endif

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

static int gillespieSSA(int *x, double *w, double *q, const int *x0, const int *checkpoints, const int rank, const int size, double *timings)
{
    // Initialize
    double t = 0.0;              // Time
    double tau;                  // Time step
    int t_index = 0;             // Time index
    double wtime;                // Wall time
    double wtime0 = MPI_Wtime(); // Wall time at start of loop
    double a0;                   // Total propensity
    int reaction;                // Reaction index

    memcpy(x, x0, Q * sizeof(int)); // Copy initial state to state vector

    // Main loop
    while (t < T)
    {
        // Check if time is larger than next checkpoint
        if (t > checkpoints[t_index])
        {
            // Write runtime to array
            wtime = MPI_Wtime() - wtime0;
            timings[t_index] += wtime;
            t_index++;
        }
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

    // Write final runtime to array
    wtime = MPI_Wtime() - wtime0;
    timings[t_index] += wtime;

    return x[0];
}

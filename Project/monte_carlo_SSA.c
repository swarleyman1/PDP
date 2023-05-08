#include <stdio.h>
#include <stdlib.h>
#include <openmpi.h>
#include <math.h>
#include <time.h>

/****************************************************************************
 * Code for the project in the course PDP at Uppsala University 2023        *
 *                                                                          *
 * An implementation of the Monte Carlo SSA with OpenMPI multiprocessing    *
 * Author: Jacob Malmenstedt                                                *
 *                                                                          *
 ****************************************************************************/

#define R 15  // Number of reactions
#define Q 7   // Number of quantities
#define T 100 // Maximum time

// Function declarations
double randdouble(double min, double max);
void update_state(int *x, int reaction);
void prop(int *x, double *w);
int select_reaction(double *q, double a0);
int gillespieSSA();

// Main function
int main(int argc, char *argv[])
{

    // Check if correct number of arguments
    if (argc != 3)
    {
        printf("Usage: %s <number of iterations> <output file>\n", argv[0]);
        return 1;
    }

    // Read command line arguments
    int n = atoi(argv[1]);    // Number of iterations per process
    char *filename = argv[2]; // Output file name

    // Initialize MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request send_request, recv_request;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Rank of current process

    // Initialize constant
    const int N = size * n; // Total number of iterations

    int local_X = (int *)malloc(n * sizeof(int)); // Local output vector

    // Start timer
    double time = MPI_Wtime();

    // Run experiments n times on each process
    for (int i = 0; i < n; i++)
    {
        local_X[i] = gillespieSSA();
    }

    // Find max and min element in local output vector
    int local_max, local_min;
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
    MPI_Reduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

    MPI_Bcast(&global_max, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global_min, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Create 20 histogram bins
    int *bins = (int *)malloc(20 * sizeof(int));
    int bin_size = (global_max - global_min) / 20;
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
    }

    // Find global bin counts
    if (rank == 0)
    {
        int *global_bin_counts = (int *)calloc(20, sizeof(int));
        MPI_Reduce(local_bin_counts, global_bin_counts, 20, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Stop timer
    time = MPI_Wtime() - time;

    // Print max time
    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("Max time: %fs\n", max_time);
    }

    // Write to file
    if (rank == 0)
    {
        FILE *fp;
        fp = fopen(filename, "w");
        fprintf(fp, "Number of iterations: %d\n", N);
        fprintf(fp, "Max time: %fs\n \n", max_time);
        fprintf(fp, "Bin \t\t Count\n");
        fprintf(fp, "------------------\n");
        for (int i = 0; i < 19; i++)
        {
            fprintf(fp, "%d - %d \t %d\n", bins[i], bins[i + 1], global_bin_counts[i]);
        }
        fprintf(fp, "%d - %d \t %d\n", bins[19], global_max, global_bin_counts[19]);
        fclose(fp);
        free(global_bin_counts);
    }

    // Free memory
    free(local_X);
    free(bins);
    free(local_bin_counts);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}

// Function definitions
double randdouble(double min, double max)
{
    // Generate random double between min and max
    double r = (double)rand() / (double)RAND_MAX;
    return min + r * (max - min);
}

// State change matrix (15x7)
void update_state(int *x, int reaction)
{
    switch (reaction)
    {
    case 0:
        x[0] += 1;
        break;
    case 1:
        x[0] -= 1;
        break;
    case 2:
        x[0] -= 1;
        x[2] += 1;
        break;
    case 3:
        x[1] += 1;
        break;
    case 4:
        x[1] -= 1;
        break;
    case 5:
        x[1] -= 1;
        x[3] += 1;
        break;
    case 6:
        x[2] -= 1;
        break;
    case 7:
        x[2] -= 1;
        x[4] += 1;
        break;
    case 8:
        x[3] -= 1;
        break;
    case 9:
        x[3] -= 1;
        x[5] += 1;
        break;
    case 10:
        x[4] -= 1;
        break;
    case 11:
        x[4] -= 1;
        x[6] += 1;
        break;
    case 12:
        x[5] -= 1;
        break;
    case 13:
        x[0] += 1;
        x[6] -= 1;
        break;
    case 14:
        x[6] -= 1;
        break;
    default:
        perror("Invalid reaction index");
        break;
    }
}

int select_reaction(double *q, double a0, int length)
{
    // Generate random double between 0 and 1
    double u2 = randdouble(0, 1);

    // Find index of first element in q that is larger than u2*a0
    int i = 0, u2a2 = u2 * a0;
    while (q[i] < u2a2 && i < length)
    {
        i++;
    }

    return i;
}

int gillespieSSA()
{
    double t = 0.0;  // Time
    double tau;      // Time step
    double a0 = 0.0; // Total propensity
    int reaction;    // Reaction index

    int *x = (int *)malloc(Q * sizeof(int));          // State vector
    x = [ 900, 900, 30, 330, 50, 270, 20 ];           // Initial state
    double *w = (double *)malloc(R * sizeof(double)); // Propensity vector
    double *q = (double *)malloc(R * sizeof(double)); // Cumulative propensity vector

    // Main loop
    while (t < T)
    {
        // Calculate propensities
        prop(x, w);

        // Calculate total propensity
        a0 = 0;
        for (int i = 0; i < R; i++)
        {
            a0 += w[i];
        }

        // Calculate time step
        tau = (1 / a0) * log(1 / randdouble(0, 1));

        // Calculate cumulative propensities
        q[0] = w[0];
        for (int i = 1; i < R; i++)
        {
            q[i] = q[i - 1] + w[i];
        }

        // Select reaction
        reaction = select_reaction(q, a0, R);

        // Update state
        update_state(x, reaction);

        // Update time
        t += tau;
    }

    return x[0]
}

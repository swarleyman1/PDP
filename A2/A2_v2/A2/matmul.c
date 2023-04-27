#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#define PRODUCE_OUTPUT_FILE

// Run with: mpirun -np 4 ./matmul input_file output_file
// number of processes must be a perfect square

/* TODO: 
Swap A and B to get correct output
Improve matrix multiplication inner loop
Remove mallocs and frees?

*/


int write_output(char *file_name, const float *output, int n, int size)
{
    // Function that writes the output in an output file
    // (transposes the matrix to get the correct output)
    FILE *file;
    if (NULL == (file = fopen(file_name, "w")))
    {
        perror("Couldn't open output file");
        return -1;
    }
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            if (0 > fprintf(file, "%.4f ", output[j + (i * n)]))
            {
                perror("Couldn't write to output file");
            }
        }
    }
    if (0 > fprintf(file, "\n"))
    {
        perror("Couldn't write to output file");
    }
    if (0 != fclose(file))
    {
        perror("Warning: couldn't close output file");
    }
    return 0;
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Status status;
    MPI_Request send_request, recv_request;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

    if (argc != 3)
    {
        printf("Usage: stencil input_file output_file\n");
        return 1;
    }

    // Initialize variables
    char *input = argv[1];
    char *output = argv[2];
    float *A, *B, *C, *A_chunk, *B_chunk, *C_chunk;
    int n;

    // Read input file
    if (rank == 0)
    {
        FILE *file;
        int x;
        file = fopen(input, "r");
        x = fscanf(file, "%d", &n);

        // Allocate memory for global matrices A, B and C
        A = malloc(n * n * sizeof(float));
        B = malloc(n * n * sizeof(float));
        C = malloc(n * n * sizeof(float));

        for (int i = 0; i < (n * n); i++)
        {
            x = fscanf(file, "%f", &A[i]);
        }
        for (int i = 0; i < (n * n); i++)
        {
            x = fscanf(file, "%f", &B[i]);
        }
    }

    // Start timer
    double start = MPI_Wtime();

    // Broadcast n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int num_rows = n / size;

    // Allocate memory for local matrices A_chunk, B_chunk and C_chunk
    A_chunk = malloc(n * num_rows * sizeof(float));
    B_chunk = malloc(n * num_rows * sizeof(float));
    C_chunk = malloc(num_rows * num_rows * sizeof(float));

    // Scatter rows of B to  the different processes
    MPI_Scatter(B, n * num_rows, MPI_FLOAT, B_chunk, n * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    for (int l = 0; l < size; l++) // Loop over chunks of A
    {
        for (int i = 0; i < n; ++i) // Create a chunk of A
        {
            for (int j = 0; j < num_rows; ++j)
            {

                if (rank == 0)
                {
                    A_chunk[j * n + i] = A[(i * n + j) + l * num_rows];
                }
            }
        }

        // broadcast the chunk of columns of A to all processes
        MPI_Bcast(A_chunk, n * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // do matrix multiplication
        memset(C_chunk, 0, num_rows * num_rows * sizeof(float));

        for (int j = 0; j < num_rows; ++j)
        {
            for (int k = 0; k < num_rows; ++k)
            {
                for (int i = 0; i < n; ++i)
                {
                    C_chunk[(j * num_rows + k)] += A_chunk[j * n + i] * B_chunk[n * k + i];
                }
            }
        }

        // gather columns of C using send and receive (not ideal but works)
        for (int q = 0; q < num_rows; q++)
        {
            if (rank == 0)
            {
                for (int r = 1; r < size; r++) // Receive from all processes (note: will not execute if serial)
                {
                    float *temp1 = malloc(num_rows * sizeof(float)); // ----------------------------- Not ideal
                    MPI_Recv(temp1, num_rows, MPI_FLOAT, r, 0, MPI_COMM_WORLD, &status);
                    for (int s = 0; s < num_rows; s++)
                    {
                        C[(r * num_rows + q * n) + (l * n * num_rows) + s] = temp1[s];
                    }
                    free(temp1);
                    for (int s = 0; s < num_rows; s++)
                    {
                        C[(q * n) + (l * n * num_rows) + s] = C_chunk[s + q * num_rows];
                    }
                }
            }
            else
            {
                float *temp2 = malloc(num_rows * sizeof(float)); // ----------------------------- Not ideal
                for (int s = 0; s < num_rows; s++)
                {
                    temp2[s] = C_chunk[q * num_rows + s];
                }
                MPI_Send(temp2, num_rows, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
                free(temp2);
            }
        }
    }

    // Stop timer
    double my_execution_time = MPI_Wtime() - start;
    double max_execution_time;

    // Find max execution time
    MPI_Reduce(&my_execution_time, &max_execution_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    { // Only process with rank 0 write output file
        printf("%lf", max_execution_time);

#ifdef PRODUCE_OUTPUT_FILE
        if (size == 1)
        {
            if (0 != write_output(output, C_chunk, n, size))
            {
                return 2;
            }
        }
        else
        {
            if (0 != write_output(output, C, n, size))
            {
                return 2;
            }
        }
#endif
    }

    // clean up
    free(A_chunk);
    free(B_chunk);
    free(C_chunk);
    if (rank == 0)
    {
        free(A);
        free(B);
        free(C);
    }
    MPI_Finalize();

    return 0;
}
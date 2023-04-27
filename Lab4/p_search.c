#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Not quite correct yet

// Run with: mpirun -np 4 ./p_search input_file search_integer
// compile with: mpicc -O3 -o p_search p_search.c

// Parallel search using MPI
// The code reads a file of integers and searches for the first occurrence of a given integer.

// Algorithm:
// 1. Read the file of integers on process 0
// 2. Broadcast the number of integers to all processes
// 3. Scatter the integers to all processes
// 4. Search for the first occurrence of the given integer in the local array
// 5. Gather the results from all processes to process 0
// 6. Print the result

int main(int argc, char *argv[])
{
    int rank, size, n, *A, *A_chunk, *result, *result_chunk, search, found = 0, found_rank = -1, found_index = -1;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Read the file of integers on process 0
    if (rank == 0)
    {
        FILE *file;
        int x;
        file = fopen(argv[1], "r");
        x = fscanf(file, "%d", &n);

        // Allocate memory for global array A
        A = malloc(n * sizeof(int));

        for (int i = 0; i < n; i++)
        {
            x = fscanf(file, "%d", &A[i]);
        }
        printf("Process %d read %d integers\n", rank, n);
        for (int i = 0; i < n; i++)
        {
            printf("%d ", A[i]);
        }
    }

    // Broadcast n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate memory for local array A_chunk
    A_chunk = malloc(n / size * sizeof(int));

    // Scatter rows of A to  the different processes
    MPI_Scatter(A, n / size, MPI_INT, A_chunk, n / size, MPI_INT, 0, MPI_COMM_WORLD);

    printf("Process %d received %d integers\n", rank, n / size);

    // Read the search integer on process 0
    if (rank == 0)
    {
        search = atoi(argv[2]);
    }

    // Broadcast search to all processes
    MPI_Bcast(&search, 1, MPI_INT, 0, MPI_COMM_WORLD);

    printf("Process %d searching for %d\n", rank, search);

    // Search for the first occurrence of the given integer in the local array
    for (int i = 0; i < n / size; i++)
    {
        if (A_chunk[i] == search)
        {
            found = 1;
            found_rank = rank;
            found_index = i;
            break;
        }
    }

    // Allocate memory for local array result_chunk
    result_chunk = malloc(3 * sizeof(int));

    // Gather the results from all processes to process 0
    MPI_Gather(&found, 1, MPI_INT, result_chunk, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int index = -1;
        if (rank == found_rank)
        {
            index = found_index + rank * n / size;
            MPI_Bcast(&index, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
    // Print the result
    if (rank == 0)
    {
        for (int i = 0; i < size; i++)
        {
            if (result_chunk[i] == 1)
            {
                found = 1;
                found_rank = i;
                // found_index = ;
                break;
            }
        }
        if (found == 1)
        {
            printf("Found %d at index %d on process %d\n", search, index, found_rank);
        }
        else
        {
            printf("Did not find %d\n", search);
        }
    }

    // Finalize MPI
    MPI_Finalize();

    // Free memory
    free(A_chunk);
    free(result_chunk);
    if (rank == 0)
    {
        free(A);
    }

    return 0;
}

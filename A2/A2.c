#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{

    if (argc != 3)
    {
        printf("Expected: Input filename and Output filename \n");
        return 0;
    }

    char *input = argv[1];
    char *output = argv[2];
    float **A, **B, **C, **A_local, **B_local;
    int n;

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int p = (int)sqrt(size);

    // Read input file
    if (rank == 0)
    {
        FILE *file;
        int x;
        file = fopen(input, "r");
        x = fscanf(file, "%d", &n);

        printf("n: %d \n", n);

        // Allocate memory for global matrices A and B
        float *A_data = malloc(n * n * sizeof(float));
        float *B_data = malloc(n * n * sizeof(float));
        A = malloc(n * sizeof(float *));
        B = malloc(n * sizeof(float *));
        for (int i = 0; i < n; i++)
        {
            A[i] = &A_data[i * n];
            B[i] = &B_data[i * n];
        }

        /*
        A = malloc(n * sizeof(float *));
        B = malloc(n * sizeof(float *));

        for (int i = 0; i < n; i++)
        {
            A[i] = malloc(n * sizeof(float));
            B[i] = malloc(n * sizeof(float));
        }
        */
       printf("scanning A \n");
        for (int i = 0; i < (n); i++)
        {
            for (int j = 0; j < (n); j++)
            {
                x = fscanf(file, "%f", &A[i][j]);
            }
        }
        printf("scanning B \n");
        for (int i = 0; i < (n); i++)
        {
            for (int j = 0; j < (n); j++)
            {
                x = fscanf(file, "%f", &B[i][j]);
            }
        }

        printf("A: \n");
        for (int i = 0; i < (n); i++)
        {
            for (int j = 0; j < (n); j++)
            {
                printf("%f \t", A[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        printf("\n B: \n");
        for (int i = 0; i < (n); i++)
        {
            for (int j = 0; j < (n); j++)
            {
                printf("%f \t", B[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /*
    A_local = malloc((n) / p * sizeof(float *));
    B_local = malloc((n) / p * sizeof(float *));

    for (int i = 0; i < n / p; i++)
    {
        A_local[i] = malloc((n) / p * sizeof(float));
        B_local[i] = malloc((n) / p * sizeof(float));
    }
    */

    // Allocate memory for local matrices A_local and B_local
    int block_size = n / p;
    float *A_local_data = malloc(block_size * block_size * sizeof(float));
    float *B_local_data = malloc(block_size * block_size * sizeof(float));
    A_local = malloc(block_size * sizeof(float *));
    B_local = malloc(block_size * sizeof(float *));
    for (int i = 0; i < block_size; i++)
    {
        A_local[i] = &A_local_data[i * block_size];
        B_local[i] = &B_local_data[i * block_size];
    }

    MPI_Datatype block_type; // Defines a block of a larger matrix
    int sizes[2] = {n, n};
    // int sizes_2[2] = {n / size, n / size};
    // int start[2] = {n / size * rank, n / size * rank};
    int start[2] = {rank % p * n / p, rank / p * n / p};
    printf("rank: %d, start: %d, %d \n", rank, start[0], start[1]);
    int sizes_2[2] = {n / p, n / p};
    MPI_Type_create_subarray(2, sizes, sizes_2, start, MPI_ORDER_C, MPI_FLOAT, &block_type);
    MPI_Type_commit(&block_type);

    MPI_Barrier(MPI_COMM_WORLD);

    // Scatter blocks of A and B to all processes
    //MPI_Scatter(A, 1, block_type, A_local, (n * n) / size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(B, 1, block_type, B_local, (n * n) / size, MPI_FLOAT, 0, MPI_COMM_WORLD);


        if (rank == 0)
        {
            for (int i = 1; i < size; i++)
            {
                // MPI_Send syntax: MPI_Send(void* data, int count, MPI_Datatype datatype, int destination, int tag, MPI_Comm communicator)
                MPI_Send(A, 1, block_type, i, 0, MPI_COMM_WORLD);
                MPI_Send(B, 1, block_type, i, 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Recv(A_local, (n * n) / size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(B_local, (n * n) / size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("A_local: \n");
        for (int i = 0; i < sqrt((n * n) / size); i++)
        {
            for (int j = 0; j < sqrt((n * n) / size); j++)
            {
                printf("%f \t", A_local[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        printf("\n B_local: \n");
        for (int i = 0; i < sqrt((n * n) / size); i++)
        {
            for (int j = 0; j < sqrt((n * n) / size); j++)
            {
                printf("%f \t", B_local[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    if (rank == 0)
    {
        //    free(A);
        //    free(B);
    }

    // free(C);

    MPI_Finalize();
}

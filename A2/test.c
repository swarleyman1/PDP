#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

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

	char *input = argv[1];
	char *output = argv[2];
	float *A, *B, *C, *A_chunk, *B_chunk, *C_chunk;
	int n;

	// int p = (int)sqrt(size);

	// Read input file
	if (rank == 0)
	{
		FILE *file;
		int x;
		file = fopen(input, "r");
		x = fscanf(file, "%d", &n);

		printf("n: %d \n", n);

		// Allocate memory for global matrices A and B

		A = malloc(n * n * sizeof(float));
		B = malloc(n * n * sizeof(float));
		C = malloc(n * n * sizeof(float));

		printf("scanning A \n");
		for (int i = 0; i < (n * n); i++)
		{
			x = fscanf(file, "%f", &A[i]);
		}
		printf("scanning B \n");
		for (int i = 0; i < (n * n); i++)
		{
			x = fscanf(file, "%f", &B[i]);
		}

		printf("A: \n");
		for (int i = 0; i < (n); i++)
		{
			for (int j = 0; j < (n); j++)
			{
				printf("%f \t", A[i * n + j]);
			}
			printf("\n");
		}
		printf("\n");
		printf("\n B: \n");
		for (int i = 0; i < (n); i++)
		{
			for (int j = 0; j < (n); j++)
			{
				printf("%f \t", B[i * n + j]);
			}
			printf("\n");
		}
		printf("\n");
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int num_rows = n / size;

	// Allocate memory for local matrices A_chunk, B_chunk and C_chunk
	A_chunk = malloc(n * num_rows * sizeof(float));
	B_chunk = malloc(n * num_rows * sizeof(float));
	C_chunk = malloc(num_rows * num_rows * sizeof(float));

	// Define the datatype for a column
	MPI_Datatype column_type;
	MPI_Type_vector(
		n,		   // count
		1,		   // blocklength
		n,		   // stride
		MPI_FLOAT, // oldtype
		&column_type);
	MPI_Type_commit(&column_type);

	// Scatter rows of b to  the different processes
	// MPI_Scatter syntax: MPI_Scatter(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm communicator)
	MPI_Scatter(B, n * num_rows , MPI_FLOAT, B_chunk, n * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);


	for (int i = 0; i < n; ++i) // are loops correct?
	{
		for (int j = 0; j < num_rows; ++j)
		{

			if (rank == 0)
			{
				A_chunk[j * n + i] = A[i * n + j];
			}
		}

		// broadcast columns of A
		// MPI_Bcast syntax: MPI_Bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator)
		MPI_Bcast(A_chunk, n * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

		// do matrix multiplication (unfinished)
		for (int j = 0; j < num_rows; ++j)
		{
			for (int k = 0; k < num_rows; ++k)
			{
				C_chunk[j * num_rows + k] += A_chunk[j * n + i] * B_chunk[i * num_rows + k]; //Wrong?
			}
		}

		// gather columns of C using send and receive
		//if(rank == 0){
		//	MPI_Send(C_chunk, n * num_rows, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
		//}else{
		//	MPI_Recv(C
		//}

		//Unnecessary
		printf("C_chunk rank %d: ", rank);
		printf("%f ", C_chunk[0]);
		printf("\n");

		C_chunk[0] = 0;


	}



	// print A_chunk and B_chunk for debugging
	printf("A_chunk rank %d: ", rank);
	for (int j = 0; j < num_rows; ++j)
	{
		for (int k = 0; k < n; ++k)
		{
			printf("%f ", A_chunk[j * n + k]);
		}
		printf("\n");
	}
	printf("B_chunk rank %d: ", rank);
	for (int j = 0; j < num_rows; ++j)
	{
		for (int k = 0; k < n; ++k)
		{
			printf("%f ", B_chunk[j * n + k]);
		}
		printf("\n");
	}

	printf("C_chunk rank %d: ", rank);
	for (int j = 0; j < num_rows; ++j)
	{
		for (int k = 0; k < num_rows; ++k)
		{
			printf("%f ", C_chunk[j * n + k]);
		}
		printf("\n");
	}

	// clean up
	free(A_chunk);
	free(B_chunk);
	if (rank == 0)
	{
		free(A);
		free(B);
		free(C);
	}
	MPI_Finalize();

	return 0;
}

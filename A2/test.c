#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int read_input(const char *file_name, float **values)
{
	FILE *file;
	if (NULL == (file = fopen(file_name, "r")))
	{
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values))
	{
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(float))))
	{
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i = 0; i < num_values; i++)
	{
		if (EOF == fscanf(file, "%f", &((*values)[i])))
		{
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file))
	{
		perror("Warning: couldn't close input file");
	}
	return num_values;
}

int main(int argc, char **argv)
{
	int rank, size;
	MPI_Status status;
	MPI_Request send_request, recv_request;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

	int n;
	int num_values;
	float *input;
	if (rank == 0)
	{
		if (argc != 3)
		{
			printf("Usage: stencil input_file output_file\n");
			return 1;
		}
		char *input_name = argv[1];
		char *output_name = argv[2];

		// Read input file
		
		
		if (0 > (num_values = read_input(input_name, &input)))
		{
			return 1;
		}

		
	}

	MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);
	n = num_values;
	float B_chunk[n * (n / size)];
	float A[n * n], B[n * n], C[n * n];
	if (rank == 0)
	{

		for (int i = 0; i < n * n; ++i)
		{
			A[i] = input[i + 1];
			B[i] = input[i + 1 + n * n];
		}

		
	}
	// Define the datatype for a column
		MPI_Datatype column_type;
		MPI_Type_vector(
			n,		   // count
			1,		   // blocklength
			n,		   // stride
			MPI_FLOAT, // oldtype
			&column_type);
		MPI_Type_commit(&column_type);

	// Send columns of b
	MPI_Scatter(B, n / size, column_type, B_chunk, n * (n / size), MPI_FLOAT, 0, MPI_COMM_WORLD);

	// Send out rows for A (unfinished)
	/*
	int rows = 1;
	float A_own[n*rows];
	if (rank == 0){
	for (int row = 0; row < n; row += rows){
		  for (int col = 0; col < n; ++col){
			A_own[col] = A[row*n+col];
		  }
	}
	}
	*/

	int num_rows = n / size;
	
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < num_rows; ++j)
		{
			float A_chunk[n * num_rows];
			if (rank == 0)
			{
				A_chunk[j * n + i] = A[i * n + j];
			}
			// broadcast rows of A
			//MPI_Bcast syntax: MPI_Bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator)
			MPI_Bcast(A_chunk, n * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

			// print A_chunk and B_chunk
			printf("A_chunk rank %d: ", rank);
			for (int j = 0; j < num_rows; ++j){
				for (int k = 0; k < n; ++k)
				{
					printf("%f ", A_chunk[j * n + k]);
				}
				printf("\n");
			}
			printf("B_chunk rank %d: ", rank);
			for (int j = 0; j < n; ++j){
				for (int k = 0; k < n; ++k)
				{
					printf("%f ", B_chunk[j * n + k]);
				}
				printf("\n");
			}

			// do matrix multiplication


			
			
		}

		// do matrix multiplication

		// gather results
	}

	// MPI_Bcast(A_own, n*rows,
	// MPI_Scatter(A, n*(n/size), MPI_FLOAT, A_own, n*(n/size), MPI_FLOAT, 0, MPI_COMM_WORLD);

	// Do matrix multiplication
	// for


	//free(B);
	return 0;
}

